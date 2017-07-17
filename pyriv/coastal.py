import networkx as nx
import geopandas as gpd
from shapely.ops import polygonize
from shapely.geometry import MultiPolygon, LineString, Point
from units import length_in_display_units
from multiprocessing import Pool
from functools import partial

#%% For testing
coastfn = '/Users/jkibele/Documents/SASAP/sasap-size-declines/RiverDistance/data/test_data/CoastLine.shp'
rivfn = '/Users/jkibele/Documents/SASAP/sasap-size-declines/RiverDistance/data/test_data/Rivers.shp'
testlinefn = '/Users/jkibele/Documents/SASAP/sasap-size-declines/RiverDistance/data/test_data/test_lines.shp'

#%% Make land polygon from coastline. Coasline is actually many linestrings.
def polygonize_coastline(coast):
    """Turn the coastline into a multipolygon.
    
    Take a geodataframe or a path to a shapefile and return a land multipolygon.
    """
    if coast.__class__.__name__ != 'GeoDataFrame':
        coast = gpd.read_file(coast)
        
    return MultiPolygon(polygonize(coast.geometry))

def coords_from_coastline(coast):
    """Extract coordinates from a coastline shapefile or geodataframe.
    
    This serves pretty much the same purpose as `extract_poly_exterior_coords`.
    We may not need to keep both versions.
    """
    if coast.__class__.__name__ != 'GeoDataFrame':
        coast = gpd.read_file(coast)
        
    clst = coast.geometry.apply(lambda l: l.coords).tolist()
    fltlst = [item for sublist in clst for item in sublist]
    return fltlst

def extract_poly_exterior_coords(geom):
    """Extract the exterior coordinates from a polygon or multipolygon.
    
    Returns
    -------
    list of tuples
        All the exterior coordinates of geom.
    
    Notes
    -----
    This code was adapted from an answer on stackexchange:
    https://gis.stackexchange.com/questions/119453/count-the-number-of-points-in-a-multipolygon-in-shapely
    """
    if geom.type == 'Polygon':
        exterior_coords = geom.exterior.coords[:]
    elif geom.type == 'MultiPolygon':
        exterior_coords = []
        for part in geom:
            epc = extract_poly_exterior_coords(part)  # Recursive call
            exterior_coords += epc[:]
    else:
        raise ValueError('Unhandled geometry type: ' + repr(geom.type))
    return exterior_coords

def graph_from_coast(coast):
    """Build a networkx graph for coastal swimmable distance calculations.
    
    """
    if coast.__class__.__name__ != 'GeoDataFrame':
        coast = gpd.read_file(coast)

def ocean_edges_for_node(node, land_obj, node_list):
    """
    multiprocessing can't work with functions that are not a top level, so I'm
    moving this out here to see if I can make it work.
    """
    ocean_edges = tuple()
    for n in node_list:
        if node <> n:
            line = LineString((node,n))
            if not land_obj.line_crosses(line):
                ocean_edges += ((node, n, {'distance': length_in_display_units(line.length)}),)
                print '.',
    return ocean_edges


#%% Land Obj
class Land(object):
    def __init__(self, coast, shrink=1.0):
        
        if coast.__class__.__name__ != 'GeoDataFrame':
            coast = gpd.read_file(coast)
            
        self.gdf = coast
        self.land_poly = polygonize_coastline(self.gdf)
        self.land_shrunk = self.land_poly.buffer(-1*shrink)
        self.coords = coords_from_coastline(self.gdf)
        self.cached_graph = None
        
    def graph(self, dump_cached=False):
        if dump_cached or not self.cached_graph:
            G = self.fresh_graph()
            G = self._add_ocean_edges_complete(G, verbose=True)
            self.cached_graph = G
        return self.cached_graph
    
    def fresh_graph(self):
        G = nx.Graph()
        G.add_nodes_from(self.coords)
        return G
    
    def line_crosses(self, line):
        """Determine whether or not the given line crosses shrunken land.
        
        """
        return line.intersects(self.land_shrunk)
    
    def _ocean_edges_for_node(self, graph, node):
        ocean_edges = tuple()
        for n in graph.nodes_iter():
            if node <> n:
                line = LineString((node,n))
                if not self.line_crosses(line):
                    ocean_edges += ((node, n, {'distance': length_in_display_units(line.length)}),)
                    print '.',
        return ocean_edges
    
    def _add_ocean_edges_complete(self, graph, n_jobs=6, verbose=False):
        if verbose:
            cnt = 1
            import time
            t0 = time.time()
            print "Starting at %s to add edges for %i nodes." % (time.asctime(time.localtime(t0)), graph.number_of_nodes() )
            edge_possibilities = graph.number_of_nodes() * (graph.number_of_nodes() -1)
            print "We'll have to look at somewhere around %i edge possibilities." % ( edge_possibilities )
            print "Node: ",
        oe_node = partial(ocean_edges_for_node, land_obj=self, node_list=graph.nodes())
        pool = Pool(processes=n_jobs)
#        o_eds = lambda n: ocean_edges_for_node(self, graph, n)
        ocean_edges = pool.map(oe_node, graph.nodes_iter())
        
#        ocean_edges = tuple()
#        for node in graph.nodes_iter():
#            if verbose:
#                print str(cnt) + ' ',
#                cnt += 1
#            ocean_edges += self._ocean_edges_for_node(graph, node)
        
        graph.add_edges_from(ocean_edges)
        if verbose:
            print "It took %i minutes to load %i edges." % ((time.time() - t0)/60, graph.number_of_edges() )
        return graph