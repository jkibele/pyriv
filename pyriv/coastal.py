import numpy as np
import networkx as nx
import geopandas as gpd
from shapely.ops import polygonize
from shapely.geometry import MultiPolygon, LineString, Point
from units import length_in_display_units
from river_graph import point_to_tuple
from multiprocessing import Pool
from functools import partial
from itertools import chain
from ast import literal_eval
import os

#%% For testing
test_data_dir = '/home/jkibele/sasap-size-declines/RiverDistance/data/test_data/'
#test_data_dir = '/Users/jkibele/Documents/SASAP/sasap-size-declines/RiverDistance/data/test_data/'
fullpath = lambda fn: os.path.join(test_data_dir, fn)
coastfn = fullpath('CoastLine.shp')
rivfn = fullpath('Rivers.shp')
testlinefn = fullpath('test_lines.shp')
rivmouthsfn = fullpath('RiverMouths.shp')
uprivpoints = fullpath('UpriverTestPoints.shp')
islandgraphmlfn = fullpath('island.graphml')
islandgpickle = fullpath('island.gpickle')


#%% Make land polygon from coastline. Coasline is actually many linestrings.
def polygonize_coastline(coast):
    """Turn the coastline into a multipolygon.
    
    Take a geodataframe or a path to a shapefile and return a land multipolygon.
    """
    if coast.__class__.__name__ != 'GeoDataFrame':
        coast = gpd.read_file(coast)
    plist = [p for p in polygonize(coast.geometry)]
    return MultiPolygon(plist)

def coords_from_coastline(coast):
    """Extract coordinates from a coastline shapefile or geodataframe.
    
    This serves pretty much the same purpose as `extract_poly_exterior_coords`.
    We may not need to keep both versions.
    """
    if coast.__class__.__name__ != 'GeoDataFrame':
        coast = gpd.read_file(coast)
    
#    print "geom type: {}".format(coast.geometry.iloc[0].type)
    if coast.geometry.iloc[0].type.lower().find('poly') > -1:
        fltlst = extract_poly_exterior_coords(coast.unary_union)
    else: # assume it's a line feature
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

def radius_filter(node, node_list, radius=5000):
    node = Point(node)
    pnts = [Point(n) for n in node_list]
    keep = [tuple(np.array(p)) for p in pnts if node.distance(p) <= radius]
    return keep

def ocean_edges_for_node(node, land_obj, node_list, radius=None):
    """
    multiprocessing can't work with functions that are not a top level, so I'm
    moving this out here to see if I can make it work.
    """
    if radius:
        node_list = radius_filter(node, node_list, radius=radius)

    ocean_edges = tuple()
    for n in node_list:
        if node <> n:
            line = LineString((node,n))
            if not land_obj.line_crosses(line):
                ocean_edges += ((node, n, {'distance': length_in_display_units(line.length)}),)
                # print '.',
    return ocean_edges

def closest_node(graph, pos):
    """
    Return the closest node to the given (x,y) position.

    Parameters
    ----------
      pos : tuple or list or shapely.geometry.point.Point
        Coordinates as (x, y), i.e., (lon, lat) not (lat, lon). If shapely
        point, will be converted to tuple.

    Returns
    -------
      tuple
        (x, y) coordinates of closest node.
    """

    nodes = np.array(graph.nodes())
    node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
    return tuple(nodes[node_pos])

def coastal_fish_distance(graph, pos0, pos1, weights='distance'):
    """Find the shortest fish-distance line from one position to another.
    
    Parameters
    ----------
    graph : nx.Graph
        Graph representation of 
    pos0 : tuple or list or shapely.geometry.point.Point
        Coordinates as (x, y), i.e., (lon, lat) not (lat, lon). If shapely
        point, will be converted to tuple.
    pos1 : tuple or list or shapely.geometry.point.Point
        Coordinates as (x, y), i.e., (lon, lat) not (lat, lon). If shapely
        point, will be converted to tuple.
    
    """
    node0 = closest_node(graph, pos0)
    node1 = closest_node(graph, pos1)
    if node0 == node1:
        pth = None
    else:
        pth = nx.shortest_path(graph, node0, node1, weight=weights)
    # convert to linestring and return
    return LineString(pth)

def cfd_to_pos_list(graph, pos0, pos_list, weights='distance'):
    """Find coastal fish distance from pos0 to a list of positions.
    """
    line_list = list()
    for pos in pos_list:
        pth = coastal_fish_distance(graph, pos0, pos, weights=weights)
        line_list.append(pth)
    return line_list


#%% Land Obj
class Land(object):
    def __init__(self, coast, graph=None, shrink=1.0):
        
        if coast.__class__.__name__ != 'GeoDataFrame':
            coast = gpd.read_file(coast)
            
        self.gdf = coast
        self.land_poly = polygonize_coastline(self.gdf)
        self.land_shrunk = self.land_poly.buffer(-1*shrink)
        self.coords = coords_from_coastline(self.gdf)
        if type(graph).__name__ == 'Graph':
            self.cached_graph = graph
        elif type(graph).__name__ == 'str':
            self.cached_graph = nx.read_gpickle(graph)
        else:
            self.cached_graph = None
        
    def graph(self, dump_cached=False, n_jobs=6, verbose=False):
        if dump_cached or not self.cached_graph:
            G = self.fresh_graph()
            G = self._add_ocean_edges_complete(G, verbose=verbose, n_jobs=n_jobs)
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
    
    def add_ocean_edges(self, nodes, n_jobs=6, radius=None):
        oe_node = partial(ocean_edges_for_node, land_obj=self, node_list=nodes, radius=radius)
        pool = Pool(processes=n_jobs)
        ocean_edges = pool.map(oe_node, self.cached_graph.nodes_iter())
        # map returns a list of tuples (one for all the edges of each node). We
        # need that flattened into a single iterable of all the edges.
        ocean_edges = chain.from_iterable(ocean_edges)
        self.cached_graph.add_edges_from(ocean_edges)
        return self.cached_graph
    
    def _add_ocean_edges_complete(self, graph, n_jobs=6, radius=None, verbose=False):
        if verbose:
            import time
            t0 = time.time()
            print "Starting at %s to add edges for %i nodes." % (time.asctime(time.localtime(t0)), graph.number_of_nodes() )
            edge_possibilities = graph.number_of_nodes() * (graph.number_of_nodes() -1)
            print "We'll have to look at somewhere around %i edge possibilities." % ( edge_possibilities )
            print "Node: ",
        oe_node = partial(ocean_edges_for_node, land_obj=self, node_list=graph.nodes(), radius=radius)
        pool = Pool(processes=n_jobs)
        ocean_edges = pool.map(oe_node, graph.nodes_iter())
        # map returns a list of tuples (one for all the edges of each node). We
        # need that flattened into a single iterable of all the edges.
        ocean_edges = chain.from_iterable(ocean_edges)
        graph.add_edges_from(ocean_edges)
        if verbose:
            print "It took %i minutes to load %i edges." % ((time.time() - t0)/60, graph.number_of_edges() )
        return graph