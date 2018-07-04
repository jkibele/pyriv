import geopandas as gpd
from tempfile import mkdtemp
import os
from .common import *
from .units import length_in_display_units
from multiprocessing import Pool
from functools import partial

from shapely.geometry import Polygon

def radius_filter(node, node_list, radius=5000):
    node = Point(node)
    pnts = [Point(n) for n in node_list]
    keep = [tuple(np.array(p)) for p in pnts if node.distance(p) <= radius]
    return keep

def ocean_edges_for_node(node, land_poly, node_list, radius=None):
    """
    Create a tuple of edges that don't cross land.

    Parameters
    ----------
    node 
        The node to calculate edges for.
    land_poly : shapely Polygon or Multipolygon
        The (shrunken) land geometry being used to create a CoastalGraph
    node_list

    radius : float
        The radius (in meters) for which to evaluate edge connections.

    Returns
    -------
    tuple
        Ocean edges to be added to graph.

    multiprocessing can't work with functions that are not a top level, so I'm
    moving this out here to see if I can make it work.
    """
    if radius is not None:
        node_list = radius_filter(node, node_list, radius=radius)

    ocean_edges = tuple()
    for n in node_list:
        if node <> n:
            line = LineString((node,n))
            if not line.intersects(land_poly):
                ocean_edges += ((node, n, {'distance': length_in_display_units(line.length)}),)
                # print '.',
    return ocean_edges

class CoastLine(nx.Graph):
    """
    This class can turn a LineString representation of a coastline into a 
    polygon representation and export a pyriv.land.Land object. Multipart
    lines will be exploded into single parts. Single part lines will be 
    arranged in order to generate polygon output.
    """
    def __init__(self, *args, **kwargs):
        """
        Build a CoastLine object.

        Parameters
        ----------
        data : input graph
            Data to initialize graph.  If data=None (default) an empty
            graph is created.  The data can be an edge list, or any
            NetworkX graph object.  If the corresponding optional Python
            packages are installed the data can also be a NumPy matrix
            or 2d ndarray, a SciPy sparse matrix, or a PyGraphviz graph.

        Returns
        -------
          A CoastLine object
        """
        super(CoastLine, self).__init__(*args, **kwargs)

    @classmethod
    def read_shp(cls, shp_fn):
        """
        Create a CoastLine oject from a shapefile. If given a shapefile with
        MultiLineString geometries, these will be exploded into singlepart
        LineStrings. Shapes are not simplified on inport, meaning that all 
        nodes are maintained in the graph.

        Parameters
        ----------
        shp_fn : str
            Either the absolute or relative path to the file or URL to be opened.
        bbox : tuple | GeoDataFrame or GeoSeries, default None
            Filter features by given bounding box, GeoSeries, or GeoDataFrame. 
            CRS mis-matches are resolved if given a GeoSeries or GeoDataFrame.
        **kwargs:
            Keyword args to be passed to the open or BytesCollection method in 
            the fiona library when opening the file. For more information on 
            possible keywords, type: import fiona; help(fiona.open)
        
        Returns
        -------
        CostLine
            A CoastLine instance
        """
        gdf = gpd.read_file(shp_fn)
        if (gdf.geom_type == 'MultiLineString').any():
            tmpdir = mkdtemp()
            shp_fn = os.path.join(tmpdir, "single_part.shp")
            sp_shp = explode(gdf).to_file(shp_fn)
        dig = nx.read_shp(shp_fn, simplify=False)
        return cls(dig)
    
    def connected_subgraphs(self):
        """
        Return an iterator of connected subgraphs. See networkx.connected_component_subgraphs
        documentation for additional information.
        """
        return nx.connected_component_subgraphs(self)
    
    def rings(self):
        """
        Return a list of 
        """
        rings = [list(nx.dfs_preorder_nodes(sg)) for sg in self.connected_subgraphs()]
        return rings
    
    def polygons(self):
        return [Polygon(r) for r in self.rings()]
    
    def poly_geodataframe(self):
        return gpd.GeoDataFrame({'geometry': self.polygons()})

    def land_object(self):
        """
        Export a pyriv.Land object for generating a CoastGraph.
        """
        return Land(self.poly_geodataframe())

class Land(gpd.GeoDataFrame):
    """
    This class will represent land geometry in order to generate a coastal graph.
    """
    def __init__(self, *args, **kwargs):
        """
        Build a RiverGraph object.
        
        Parameters
        ----------
        shrink : float, optional
            The distance (in meters) that land will be shrunk by for calculating
            edge intersection with land. If zero, paths that touch the coast will
            be eliminated. The default (1.0 m) allows for paths directly along the 
            coastline.
            
        Returns
        -------
          A Land object
        """
        # have to use `object.__setattr___` because GeoPandas and Pandas handle attr
        # creation differently to allow for columns to be attributes.
        object.__setattr__(self, 'shrink', kwargs.pop('shrink', 1.0))
        super(Land, self).__init__(*args, **kwargs)
        self._init_properties()

    def _init_properties(self):
        land_geom = self.unary_union
        land_shrunk = land_geom.buffer(-1 * self.shrink)
        object.__setattr__(self, 'land_geom', land_geom)
        object.__setattr__(self, 'land_shrunk', land_shrunk)

    def line_crosses(self, line):
        """
        Test if the line crosses land. We use the `land_shrunk` geometry to determine this.
        Shrinking the land allows lines that are right along the coast to count as not
        crossing land. This is typically preferable.

        Parameters
        ----------
        line : shapely LineString
            The line to test for crossing land.

        Returns
        -------
        bool
            True if the line crosses land, False if it does not.
        """
        return line.intersects(self.land_shrunk)

    def add_ocean_edges_complete(self, graph, n_jobs=6, radius=None, verbose=False):
        if verbose:
            import time
            t0 = time.time()
            print "Starting at %s to add edges for %i nodes." % (time.asctime(time.localtime(t0)), graph.number_of_nodes() )
            edge_possibilities = graph.number_of_nodes() * (graph.number_of_nodes() -1)
            print "We'll have to look at somewhere around %i edge possibilities." % ( edge_possibilities )
            print "Node: ",
        oe_node = partial(ocean_edges_for_node, land_poly=self.land_shrunk, node_list=graph.nodes(), radius=radius)
        pool = Pool(processes=n_jobs)
        ocean_edges = pool.map(oe_node, graph.nodes_iter())
        # map returns a list of tuples (one for all the edges of each node). We
        # need that flattened into a single iterable of all the edges.
        ocean_edges = chain.from_iterable(ocean_edges)
        graph.add_edges_from(ocean_edges)
        if verbose:
            print "It took %i minutes to load %i edges." % ((time.time() - t0)/60, graph.number_of_edges() )
        return graph

    
