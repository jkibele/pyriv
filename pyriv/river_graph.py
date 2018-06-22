import networkx as nx
import numpy as np
import json
from shapely.geometry import LineString, point, Point
import geopandas as gpd
from common import *

class RiverGraph(nx.DiGraph):
    """
    A graph representation of a river network.
    
    See RiverGraph.__init__ docstring
    """
    def __init__(self, *args, **kwargs):
        """
        Build a RiverGraph object.
        
        Parameters
        ----------
          data : networkx.DiGraph or path to shapefile
            This is what the river network graph will be built from. If it's a 
            path to a shapefile, it'll be converted to a networkx graph.
          coastline : path to shapefile or a geopandas.GeoDataFrame object
            This can be a line or polygon representation of the coastline
          riv_mouth_buffer : float
            A distance in map units (generally inteded to be meters). River
            deadends within this distance of the coast will be considered to
            be river mouths. The default is 1.5 meters.
            
        Returns
        -------
          A RiverGraph object
        """

        #if you have problems with typing, rememeber self.as_super 
        #magic word here
        if "coastline" in kwargs.keys():
            coastline_shp = kwargs["coastline"]
            self.coastline, self.crs = get_coastline_geom(coastline_shp)
        else:
            self.coastline = None
            self.crs = None
        self._river_mouths_cache = None
        self._inland_deadends_cache = None
        self._deadends_cache = None
        if "data" in kwargs.keys():
            if type(kwargs["data"]) == str:
                # handle a shapefile path
                kwargs["data"] = nx.read_shp(kwargs["data"])
                
        if "riv_mouth_buffer" in kwargs.keys():
            self.riv_mouth_buffer = kwargs["riv_mouth_buffer"]
        else:
            self.riv_mouth_buffer = 1.5
        
        self = super(RiverGraph, self).__init__(*args, **kwargs)

    @property
    def river_mouths(self):
        if not self._river_mouths_cache:
            ddf = self.deadend_gdf()
            self._river_mouths_cache = ddf[ddf.is_coastal].geometry.apply(lambda g: tuple(np.array(g))).tolist()
            self._inland_deadends_cache = ddf[~ddf.is_coastal].geometry.apply(lambda g: tuple(np.array(g))).tolist()
        return self._river_mouths_cache

    @property
    def inland_deadends(self):
        if not self._river_mouths_cache:
            ddf = self.deadend_gdf()
            self._river_mouths_cache = ddf[ddf.is_coastal].geometry.apply(lambda g: tuple(np.array(g))).tolist()
            self._inland_deadends_cache = ddf[~ddf.is_coastal].geometry.apply(lambda g: tuple(np.array(g))).tolist()
        return self._inland_deadends_cache
    
    @property
    def river_deadends(self):
        if not self._deadends_cache:
            self._deadends_cache = self.deadends()
        return self._deadends_cache
    
    def delete_cache(self):
        """
        Delete `_river_mouths_cache`, `_inland_deadends_cache`, and `_deadends_cache`.
        """
        self._river_mouths_cache = None
        self._inland_deadends_cache = None
        self._deadends_cache = None

    def closest_node(self, pos):
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

        nodes = np.array(self.nodes())
        node_pos = np.argmin(np.sum((nodes - pos)**2, axis=1))
        return tuple(nodes[node_pos])

    def get_path(self, n0, n1):
        """
        If n0 and n1 are adjacent connected nodes in the graph, this function
        return an array of point coordinates along the river linking
        these two nodes.
        """
        return np.array(json.loads(self[n0][n1]['Json'])['coordinates'])

    def edge_distance(self, ed):
        """
        Calculate edge distance by finding all vertices between nodes,
        creating a linestring, and return its distance. With the NHD data set
        this isn't necessary, because the edges have an attribute called
        'LengthKM'.

        Parameters
        ----------
          ed : tuple or array-like
            The coordinates of the nodes defining the edge. (node0, node1)
            a.k.a. ((x0,y0),(x1,y1)). These must be the exact coords of
            the nodes.

        Returns
        -------
          length : float
            The edge length in units that are dependent on the projection.
        """
        pth = self.get_path(*ed)
        return LineString(pth).length

    def weight_edges(self):
        """
        Calculate distance for each edge in a graph and return a copy of the
        graph with `distance` assigned to each edge. With the NHD data set this
        isn't necessary, because the edges have an attribute called 'LengthKM'.
        """
        sg = self
        for ed in sg.edges_iter():
            sg[ed[0]][ed[1]]['distance'] = self.edge_distance(ed)
        return sg

    def get_full_path(self, short_path):
        """
        This will take a path consisting of just nodes, and return the
        full path that contains all the vertices between nodes as well.
        """
        pnts = []
        for i, pnt in enumerate(short_path[:-1]):
            p = self.get_path(pnt, short_path[i+1])
            pnts.append(p)
        return np.vstack(pnts)

    def shortest_full_path(self, pos0, pos1, weights='distance'):
        """
        Find the shortest full path between 2 positions. This will find
        the nodes closest to the given positions, the nodes between those
        nodes for the shortest path, and fill in all the available
        vertices in the linestring.

        Parameters
        ----------
          pos0 : tuple or array-like
            x,y (lon,lat) coordinates for the start point
          pos1 : tuple or array-like
            x,y (lon,lat) coordinates for the end point
          sg : NetworkX Graph or DiGraph object
            The Graph

        Returns
        -------
          shapely.geometry.linestring.LineString
            A geometry representing the shortest full path.
        """
        # find the nodes closest to given (x,y) positions
        n0 = self.closest_node(pos0)
        n1 = self.closest_node(pos1)
        # find the nodes that define the shortest path between nodes
        pth = nx.shortest_path(self, n0, n1, weight=weights)
        # fill in the missing vertices
        pth = self.get_full_path(pth)
        # convert to linestring and return
        return LineString(pth)

    def is_deadend(self, node):
        if self.successors(tuple(node)):
            return False
        else:
            return True

    def deadends(self):
        """
        For a directed graph, find the nodes that have in edges but no out edges.
        """
        # get a dictionary of { node: out_degree }
        degdict = self.out_degree(self.nodes())
        # convert to a 2 column array of node and out_degree
        degarr = np.array(degdict.items())
        # get a boolean index of rows where out_degree==0
        dead_row_ind = (degarr[:,1]==0)
        # use that index to get the nodes
        dead_nodes = degarr[:,0][dead_row_ind]
        return tuple(dead_nodes)

    def deadend_gdf(self, dist=None):
        """
        Create a point geodataframe of deadend nodes attributed with wheter each 
        deadend is coastal or inland.

        Parameters
        ----------
          dist : float
            Threshold distance from the coast for a node to be considered coastal.
            Distance units are the units of the geographic projection. In the case
            of Alaska Albers the unit is meters. If `dist` is left as the default
            value `None`, the value from the `RiverGraph.riv_mouth_buffer` 
            attribute will be used (default is 1.5).

        Return
        ------
          geodataframe
            Geodataframe of point geometry attributed with boolean `is_coastal` and
            string `end_type` with values of 'Coastal' and 'Inland'.
        """
        if dist is None:
            dist = self.riv_mouth_buffer
        dnodes = self.deadends()
        ddf = gpd.GeoDataFrame({'geometry': [Point(n) for n in dnodes]})
        if self.crs:
            ddf.crs = self.crs
        if self.coastline is not None:
            ddf['is_coastal'] = ddf.distance(self.coastline) <= dist
            ddf['end_type'] = ddf.is_coastal.map({True: 'Coastal', False: 'Inland'})
        return ddf

    def is_rivermouth(self, node):
        """
        Return `True` if the given node is coastal.
        """
        node = tuple(node)
        if node in self.river_mouths:
            return True
        else:
            return False

    def reachable_nodes(self, start_node):
        """
        Return an array of reachable downstream nodes from `start_node` in
        the directed graph `sg`.
        """
        reach_set = nx.descendants(self, start_node)
        return list(reach_set)

    def reachable_subgraph(self, start_node, include_start=True):
        """
        Return a subgraph of downstream edges and nodes that can be reached by
        going downstream from `start_node`.
        """
        rns = self.reachable_nodes(start_node)
        if include_start:
            rns.append(start_node)
        return nx.subgraph(self, rns)
    
    def downstream_deadends(self, start_node):
        rn = self.reachable_nodes(start_node)
        return [n for n in rn if n in self.river_deadends]

    def downstream_rivermouths(self, start_node):
        """
        Given a start_node, return only the reachable nodes that are coastal.
        """
        rn = self.reachable_nodes(start_node)
        return [n for n in rn if n in self.river_mouths]

    def has_rivermouth(self, node_list):
        """
        Given a list of nodes, return `True` if at least one node is a river 
        mouth. Otherwise, return `False`.
        """
        return [n for n in node_list if n in self.river_mouths].any()
    
    def prune_network(self, verbose=False):
        """
        Remove subgraphs of the network that do not connect to the coastline.
        """
        # weakly_connected_component_subgraphs doesn't work with the 
        # RiverGraph subclass, so we have to switch it back to DiGraph
        sg = self.graph
        coast_n = 0
        noncoast_n = 0
        g_list = []
        cps = nx.weakly_connected_component_subgraphs(sg)
        for cp in cps:
            if has_rivermouth(cp.nodes(), sg):
                g_list.append(cp)
                coast_n += 1
            else:
                noncoast_n += 1
        if verbose:
            print "{} graphs with coastal nodes, {} without.".format(coast_n, noncoast_n)
        return RiverGraph(data=nx.compose_all(g_list), coastline_shp=self.coast_fn)

    def paths_to_coast(self, start_node, weights='distance'):
        """
        Return a list of LineStrings that represent the shortest respective
        paths to reachable `terminal_coastal_nodes`.
        """
        path_list = []
        for tcn in self.downstream_rivermouths(start_node):
            sfp = self.shortest_full_path(start_node, tcn, weights=weights)
            path_list.append(sfp)
        return path_list
    
    def paths_to_deadends(self, start_node, weights='distance'):
        """
        Return a list of LineStrings that represent the shortest respective
        paths to reachable `terminal_coastal_nodes`.
        """
        path_list = []
        for tcn in self.downstream_deadends(start_node):
            sfp = self.shortest_full_path(start_node, tcn, weights=weights)
            path_list.append(sfp)
        return path_list

    def shortest_path_to_coast(self, start_node, weights='distance'):
        """
        Return the list of nodes that constitutes the shortest path to the coast.
        """
        ddict = {}
        # in most cases, there will only be one path to the coast, but we need
        # to handle the other cases too.
        for lnst in self.paths_to_coast(start_node, weights=weights):
            # populate dictionary with path lengths as keys and LineStrings as
            # values
            ddict[lnst.length] = lnst
        # handle the case where there's no path
        if not ddict:
            # Return an empty LineString. Length will be 0
            result = LineString()
        else:
            result = ddict[min(ddict.keys())]
        return result
    
    def shortest_path_to_deadend(self, start_node, weights='distance'):
        """
        Return the list of nodes that constitutes the shortest path to the coast.
        """
        ddict = {}
        # in most cases, there will only be one path to the dead, but we need
        # to handle the other cases too.
        for lnst in self.paths_to_deadends(start_node, weights=weights):
            # populate dictionary with path lengths as keys and LineStrings as
            # values
            ddict[lnst.length] = lnst
        # handle the case where there's no path
        if not ddict:
            # Return an empty LineString. Length will be 0
            result = LineString()
        else:
            result = ddict[min(ddict.keys())]
        return result

    def river_distances(self, pnts_gdf, node_distance=True):
        """
        Assumes projection units are meters for the moment.
        """
        if pnts_gdf.__class__.__name__ != 'GeoDataFrame':
            pnts = gpd.read_file(shp_file)
        else:
            pnts = pnts_gdf
        # function to find closest node for a point geometry
        cl_nd = lambda p: self.closest_node(point_to_tuple(p))
        if self.coastline is None:
            path_find = lambda p: self.shortest_path_to_deadend(cl_nd(p))
        else:    
            path_find = lambda p: self.shortest_path_to_coast(cl_nd(p))
        pnts['path'] = pnts.geometry.apply(path_find)
        pth_dist = lambda p: p.length * 0.001
        pnts['rivdist_km'] = pnts.path.apply(pth_dist)
        if node_distance:
            pnts['nearest_node'] = pnts.geometry.apply(lambda p: Point(cl_nd(p)))
            d = {}
            for i, row in pnts.iterrows():
                d[i] = row['geometry'].distance(row['nearest_node'])
            dser = pd.Series(d, name='node_dist')
            pnts = pnts.join(dser)
        return pnts

    def plot(self, **kwargs):
        pos = dict(zip(self.nodes(), self.nodes()))
        return nx.draw_networkx_edges(self, pos, **kwargs)
