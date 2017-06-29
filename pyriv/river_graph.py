import networkx as nx
import numpy as np
import json
from shapely.geometry import LineString, point, Point
import geopandas as gpd

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
    cldf = gpd.read_file(shape_fn)
    return cldf.unary_union

class RiverGraph(nx.DiGraph):
    """
    A graph representation of a river network.
    """
    def __init__(self, coastline_shp=None, *args, **kwargs):
        """
        To make a RiverGraph from a graph, RiverGraph(data=graph)
        """
        if coastline_shp:
            self.coastline = get_coastline_geom(coastline_shp)
        else:
            self.coastline = None
        self = super(RiverGraph, self).__init__(*args, **kwargs)

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

    def shortest_full_path(self, pos0, pos1, weights='LengthKM'):
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

    def deadend_gdf(self, epsg=None, dist=1.5):
        dnodes = self.deadends()
        ddf = gpd.GeoDataFrame({'geometry': [Point(n) for n in dnodes]})
        if epsg:
            ddf.crs = {'init' :'epsg:{}'.format(epsg)}
        if self.coastline:
            ddf['is_coastal'] = ddf.distance(self.coastline) <= dist
            ddf['end_type'] = ddf.is_coastal.map({True: 'Coastal', False: 'Inland'})
        return ddf

    def is_coastal_node(self, node):
        """
        Return `True` if the given node is coastal. We call a node coastal
        if it is connected to an edge that has an `FCode` of 56600. According
        to the NHDFlowline feature type list (nhd.usgs.gov), that `FCode`
        means "coastline".
        """
        node = tuple(node)
        # get the edges with FCodes
        edges = self.edges(node, data='FCode')
        # get a list of FCodes from edges connected to the node
        fc_list = [fc for e1, e2, fc in edges]
        if self.fcode in fc_list:
            result = True
        else:
            result = False
        return result

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

    def reachable_coast_nodes(self, start_node):
        """
        Given a start_node, return only the reachable nodes that are coastal.
        """
        rn = np.array(self.reachable_nodes(start_node))
        return rn[np.apply_along_axis(self.is_coastal_node, 1, rn)].tolist()

    def has_coast_node(self, node_list):
        """
        Given a list of nodes, return `True` if at least one node is
        coastal. Otherwise, return `False`.
        """
        return np.apply_along_axis(self.is_coastal_node, 1, np.array(node_list)).any()

    def only_coastal_predecessors(self, node):
        """
        For a given node in the RiverGraph, determine if a node has one or more
        coastal node predecessors and no other predecessors.
        """
        node = tuple(node)
        preds = self.predecessors(node)
        if preds:
            only_coastal = np.array([self.is_coastal_node(nd) for nd in preds]).all()
        else:
            only_coastal = False
        return only_coastal

    def __only_coastal_index(self, node_list):
        """
        Generate a boolean index for node_list that indicates whether or not each
        node has predecessors that are exclusively coastal nodes.
        """
        nla = np.array(node_list)
        ind = np.apply_along_axis(self.only_coastal_predecessors, 1, nla)
        return ind

    def redundant_coastal_nodes(self, node_list):
        """
        Given a list of nodes, return the subset that are only reachable from
        other coastal nodes. When trying to find the shortest path to the
        coast, these nodes are not of interest.
        """
        only_coastal_bool = self.__only_coastal_index(nla)
        nla = np.array(node_list)
        return nla[only_coastal_bool].tolist()

    def terminal_coastal_nodes(self, start_node):
        """
        For a given `start_node`, return a list of nodes where the downstream
        flow meets the coastline.
        """
        rcns = self.reachable_coast_nodes(start_node)
        oci = self.__only_coastal_index(rcns)
        # Return the reachable coastal nodes that do NOT have only coastal
        # predecessors. i.e., coastal nodes with non-coastal flow into them.
        return np.array(rcns)[~oci].tolist()

    def paths_to_coast(self, start_node, weights='LengthKM'):
        """
        Return a list of LineStrings that represent the shortest respective
        paths to reachable `terminal_coastal_nodes`.
        """
        path_list = []
        # this makes it much faster
        subg = self.reachable_subgraph(start_node)
        for tcn in subg.terminal_coastal_nodes(start_node):
            sfp = subg.shortest_full_path(start_node, tcn, weights=weights)
            path_list.append(sfp)
        return path_list

    def shortest_path_to_coast(self, start_node, weights='LengthKM'):
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
            result = None
        else:
            result = ddict[min(ddict.keys())]
        return result

    def plot(self, **kwargs):
        pos = dict(zip(self.nodes(), self.nodes()))
        return nx.draw_networkx_edges(self, pos, **kwargs)
