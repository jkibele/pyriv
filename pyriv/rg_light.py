import networkx as nx
import numpy as np
import json
from shapely import ops
from shapely.geometry import LineString, point, Point, LinearRing
import geopandas as gpd
import pandas as pd

def json_linestring_reverse(ls_json):
    ed_attr_dict = json.loads(ls_json)
    ed_attr_dict["coordinates"] = ed_attr_dict["coordinates"][::-1]
    rev_json = json.dumps(ed_attr_dict)
    return rev_json

def full_reverse(G):
    """
    This will reverse the linestring path between nodes as well as the
    order of the nodes.
    """
    G = G.reverse()
    for n0,n1 in G.edges_iter():
        try:
            ls_json = G[n0][n1]['Json']
            G[n0][n1]['Json'] = json_linestring_reverse(ls_json)
        except KeyError:
            # This means there's no json path so we don't need to 
            # reverse it.
            pass
    return G

def add_reverse(G):
    """
    Make a fully reversed directional copy and add it to the original. 
    This will let us find the proper path up and down the same river.
    """
    rG = full_reverse(G)
    return nx.compose(G, rG)


def path_completion(land_df, path_geom):
    """
    Don't pass in a point that's already outside the land polygon.
    """
    if type(path_geom).__name__ == 'LineString':
        riv_mouth_pnt = Point(path_geom.coords[-1])
        ncp = nearest_coast_pnt(land_df, riv_mouth_pnt)
        ret_path = ops.linemerge((path_geom, LineString((riv_mouth_pnt, ncp))))
    else: # sometimes I want to just pass in a point
        ncp = nearest_coast_pnt(land_df, path_geom)
        ret_path = LineString((path_geom, ncp))
    return ret_path
    

def nearest_coast_pnt(land_df, riv_mouth_pnt, return_dist=False):
    """
    in pseudo-code:
        if riv mouth point is outside land poly, just return it.
        if in land poly find nearest point on outer ring of that poly
    """
    in_poly = land_df.intersects(riv_mouth_pnt)
    if not in_poly.any():
        c_pnt = riv_mouth_pnt
        c_dist = 0.0
    else:
        lglr = LinearRing(land_df[in_poly].geometry.iloc[0].exterior)
        c_dist = lglr.project(riv_mouth_pnt)
        c_pnt = lglr.interpolate(c_dist)
    
    if return_dist:
        return riv_mouth_pnt.distance(c_pnt)
    else:
        return c_pnt

def nearest_river_pnt(riv_df, pnt, return_dist=False, threshold=None):
    """
    Find nearest point on a line in riv_df.
    """
    nearest_riv = riv_df.loc[riv_df.distance(pnt).argmin()]
    rivdist = nearest_riv.geometry.distance(pnt)
    
    r_geom = nearest_riv.geometry
    r_dist = r_geom.project(pnt)
    r_pnt = r_geom.interpolate(r_dist)
    
    if return_dist:
        return pnt.distance(r_pnt)
    elif threshold:
        if pnt.distance(r_pnt) > threshold:
            return None
        else:
            return r_pnt
    else:
        return r_pnt

def deadend_distances(shp_file, riv_graph, node_distance=False):
    """
    Assumes projection units are meters for the moment.
    """
    if type(shp_file).__name__ == 'GeoDataFrame' :
         pnts = shp_file
    else:
         pnts = gpd.read_file(shp_file)
    # function to find closest node for a point geometry
    cl_nd = lambda p: riv_graph.closest_node(point_to_tuple(p))
    path_find = lambda p: riv_graph.shortest_path_to_deadend(cl_nd(p))
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

def point_to_tuple(g):
    """
    Convert `shapely.geometry.point.Point` or geopandas geoseries of length 1
    to x,y coordinate tuple. If given a tuple, list, or numpy array a tuple will
    be returned.
    """
    # this will handle geoseries
    if g.__class__.__name__ == 'GeoSeries':
        if len(g) == 1:
            g = g.iloc[0]
        else:
            raise Exception("Can not convert {} length geoseries to tuple. Only length 1.".format(len(g)))
    elif g.__class__.__name__ == 'GeoDataFrame':
        if len(g) == 1:
            g = g.geometry.iloc[0]
        else:
            raise Exception("Can not convert {} length geodataframe to tuple. Only length 1.".format(len(g)))

    if g.__class__.__name__ in ['tuple', 'list', 'ndarray']:
        result = tuple(g)
    else:
        # this will handle shapely `Point` geometries
        result = tuple(np.array(g))
    return result

def get_coastline_geom(shape_fn):
    """
    Read a shapefile, run unary_union on the geometries and return the resulting
    geometry. In the case of a coastline, this will be a multilinestring of the
    coast.

    Parameters
    ----------
      shape_fn : string
        The filepath to a line shapefile. Coastline polygons won't work.

    Returns
    -------
      shapely.geometry.MultiLinestring
        Just the geometry. Ready to use for distance calculations.
    """
    if shape_fn:
        cldf = gpd.read_file(shape_fn)
        return cldf.unary_union, cldf.crs
    else:
        return None, None

class RiverGraph(nx.DiGraph):
    """
    A graph representation of a river network.
    """
    def __init__(self, *args, **kwargs):
        """
        To make a RiverGraph from a graph, RiverGraph(data=graph)
        """

        #if you have problems with typing, rememeber self.as_super 
        #magic word here
#        self.coastline, self.crs = get_coastline_geom(coastline_shp) 
#        self._river_mouths_cache = None
#        self._inland_deadends_cache = None
        self._deadends_cache = None
        #print kwargs['data']
        #print "BLAHHNSIPUJQO"
        #self = super.DiGraph = kwargs['data']
        # --->>>>>>>> this line causes the graph to be generated twice!!!
        #               b/c of how nx works  
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

    def deadend_gdf(self, dist=1.5):
        """
        Create a point geodataframe of deadend nodes attributed with wheter each 
        deadend is coastal or inland.

        Parameters
        ----------
          dist : float
            Threshold distance from the coast for a node to be considered coastal.
            Distance units are the units of the geographic projection. In the case
            of Alaska Albers the unit is meters.

        Return
        ------
          geodataframe
            Geodataframe of point geometry attributed with boolean `is_coastal` and
            string `end_type` with values of 'Coastal' and 'Inland'.
        """
        dnodes = self.deadends()
        ddf = gpd.GeoDataFrame({'geometry': [Point(n) for n in dnodes]})
        if self.crs:
            ddf.crs = self.crs
        ddf['is_coastal'] = ddf.distance(self.coastline) <= dist
        ddf['end_type'] = ddf.is_coastal.map({True: 'Coastal', False: 'Inland'})
        return ddf

    def is_rivermouth(self, node):
        """
        Return `True` if the given node is coastal.
        """
        node = tuple(node)
        # get the edges with FCodes
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


    # I think I'll delete this. Just commenting out for now
    # def only_coastal_predecessors(self, node):
    #     """
    #     For a given node in the RiverGraph, determine if a node has one or more
    #     coastal node predecessors and no other predecessors.
    #     """
    #     node = tuple(node)
    #     preds = self.predecessors(node)
    #     if preds:
    #         only_coastal = np.array([self.is_coastal_node(nd) for nd in preds]).all()
    #     else:
    #         only_coastal = False
    #     return only_coastal

    # def __only_coastal_index(self, node_list):
    #     """
    #     Generate a boolean index for node_list that indicates whether or not each
    #     node has predecessors that are exclusively coastal nodes.
    #     """
    #     nla = np.array(node_list)
    #     ind = np.apply_along_axis(self.only_coastal_predecessors, 1, nla)
    #     return ind

    # def redundant_coastal_nodes(self, node_list):
    #     """
    #     Given a list of nodes, return the subset that are only reachable from
    #     other coastal nodes. When trying to find the shortest path to the
    #     coast, these nodes are not of interest.
    #     """
    #     only_coastal_bool = self.__only_coastal_index(nla)
    #     nla = np.array(node_list)
    #     return nla[only_coastal_bool].tolist()

    # def terminal_coastal_nodes(self, start_node):
    #     """
    #     For a given `start_node`, return a list of nodes where the downstream
    #     flow meets the coastline.
    #     """
    #     rcns = self.reachable_coast_nodes(start_node)
    #     oci = self.__only_coastal_index(rcns)
    #     # Return the reachable coastal nodes that do NOT have only coastal
    #     # predecessors. i.e., coastal nodes with non-coastal flow into them.
    #     return np.array(rcns)[~oci].tolist()

    def paths_to_coast(self, start_node, weights='distance'):
        """
        Return a list of LineStrings that represent the shortest respective
        paths to reachable `terminal_coastal_nodes`.
        """
        path_list = []
        # this makes it much faster
        # subg = self.reachable_subgraph(start_node)
        # actually, no. I don't think that makes it faster
        # ...and it doesn't work with reclass of DiGraph
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

    def plot(self, **kwargs):
        pos = dict(zip(self.nodes(), self.nodes()))
        return nx.draw_networkx_edges(self, pos, **kwargs)
