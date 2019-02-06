import networkx as nx
from networkx import NetworkXNoPath
import pandas as pd
import geopandas as gpd
import numpy as np
import json
from itertools import combinations
from shapely.geometry import LineString

def mad_matrix(mad_dataframe):
    """Turn a mad_dataframe (output of `mad_df`) into a matrix format.

    """
    mdf = mad_dataframe
    labels = mdf['from'].append(mdf['to']).unique()
    labels.sort()
    dist_rows = []
    for sc_row in labels:
        drow = []
        for sc_col in labels:
            if sc_col == sc_row:
                drow.append(0.0)
            else:
                dval_ind = ((mdf['from']==sc_row) & (mdf['to']==sc_col))
                if not dval_ind.any():
                    dval_ind = ((mdf['from']==sc_col) & (mdf['to']==sc_row))
                    if not dval_ind.any():
                        print "There's a problem with {} and {}.".format(sc_row, sc_col)
                drow.append(mdf.loc[dval_ind, 'dist_km'].item())
        dist_rows.append(drow)

    return pd.DataFrame(dist_rows, index=labels, columns=labels)


class MadGraph(nx.DiGraph):
    """
    This class represents a network graph of routes that can be traversed without
    crossing land to find Minimum Aquatic Distance (MAD).

    While it is technically a directed graph (because it subclasses DiGraph), it
    behaves as n undirected graph. When constructed, the reverse of each edge is
    included in the graph. This is necessary because, in the river portions of the
    graph, the line segments between nodes are stored as edge attributes. If we
    simply change the directed graph to undirected, the line segment geometries
    are still ordered, so traversing the edge backwards yeilds weird and wrong
    results.
    """
    def __init__(self, *args, **kwargs):
        """
        Build a MadGraph object.

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
          A MadGraph object
        """
        super(MadGraph, self).__init__(*args, **kwargs)

    def write_gpickle(self, file_path):
        nx.write_gpickle(self, file_path)

    @classmethod
    def from_gpickle(cls, file_path):
        return nx.read_gpickle(file_path)

    def round_nodes(self, decimal_places=3):
        """
        Round all node coordinates in place. Assuming meters as units, the
        default value of 3 decimal places means that nodes will be rounded
        to the nearest millimeter.
        """
        myround = lambda f: round(f, decimal_places)
        tround = lambda t: tuple([myround(f) for f in t])
        self = nx.relabel_nodes(self, tround, copy=False)
        return self

    @property
    def edge_linestrings(self):
        """
        Create a shapely LineString geometry for each edge and return a list
        of those geometries.
        """
        pths = [LineString(e) for e in self.edges()]
        return pths

    @property
    def edge_geodataframe(self):
        """
        Create a shapely LineString geometry for each edge and return a
        geopandas.GeoDataFrame containing those geometries.
        """
        pgdf = gpd.GeoDataFrame({'geometry': self.edge_linestrings})
        return pgdf

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
        node_pos = np.argmin(np.sum((nodes - pos) ** 2, axis=1))
        return tuple(nodes[node_pos])

    def get_path(self, n0, n1):
        """
        If n0 and n1 are adjacent connected nodes in the graph, this function
        return an array of point coordinates along the river linking
        these two nodes.
        """
        try:
            linepath = np.array(json.loads(self[n0][n1]['Json'])['coordinates'])
        except KeyError:
            linepath = np.array((n0, n1))
        return linepath

    def get_full_path(self, short_path):
        """
        This will take a path consisting of just nodes, and return the
        full path that contains all the vertices between nodes as well.
        """
        pnts = []
        for i, pnt in enumerate(short_path[:-1]):
            p = self.get_path(pnt, short_path[i + 1])
            pnts.append(p)
        return np.vstack(pnts)

    def shortest_full_path(self, node0, node1, weights='distance'):
        """
        Find the shortest full path between 2 positions. This will find
        the nodes closest to the given positions, the nodes between those
        nodes for the shortest path, and fill in all the available
        vertices in the linestring.

        Parameters
        ----------
          node0 : tuple or array-like
            x,y (lon,lat) coordinates for the start point
          node1 : tuple or array-like
            x,y (lon,lat) coordinates for the end point
          self : NetworkX Graph or DiGraph object
            The Graph

        Returns
        -------
          shapely.geometry.linestring.LineString
            A geometry representing the shortest full path.
        """
        # find the nodes that define the shortest path between nodes
        pth = nx.shortest_path(self, node0, node1, weight=weights)
        # fill in the missing vertices
        pth = self.get_full_path(pth)
        # convert to linestring and return
        return LineString(pth)

    def coastal_fish_distance(self, pos0, pos1, weights='distance', raise_fail=True):
        """Find the shortest fish-distance (aka minimum aquatic distance) line from
        one position to another.

        Parameters
        ----------
          self : nx.Graph
            Graph representation of the combined river and coastal network.
          pos0 : tuple or list or shapely.geometry.point.Point
            Coordinates as (x, y), i.e., (lon, lat) not (lat, lon). If shapely
            point, will be converted to tuple.
          pos1 : tuple or list or shapely.geometry.point.Point
            Coordinates as (x, y), i.e., (lon, lat) not (lat, lon). If shapely
            point, will be converted to tuple.

        """
        node0 = self.closest_node(pos0)
        node1 = self.closest_node(pos1)
        if node0 == node1:
            pth = None
        else:
            try:
                pth = self.shortest_full_path(node0, node1, weights=weights)
            except NetworkXNoPath:
                pth = None
                if raise_fail:
                    raise NetworkXNoPath("No network path between these nodes.")
        # convert to linestring and return
        return LineString(pth)

    def mad_df(self, locations, label_column, weights='distance', raise_fail=True):
        """Create a minimum aquatic distance (MAD) geodataframe for all points
        in `locations`.

        Parameters
        ----------
          self : nx.Graph
            Graph representation of the combined river and coastal network.
          locations : geopandas.GeoDataFrame or a string file path
            A point GeoDataFrame (or file path that can be opened by
            `geopandas.read_file`) representing the locations that you want
            a MAD matrix for.
          label_column : string
            The name of the column in `locations` that contains labels you
            want to use in the output.

        Returns
        -------
          geopandas.geodataframe
            A geodataframe containing MAD paths between each point in
            `locations`. In addtion to the path geometries, it will also
            have the following columns:
              from: the `label_column` value for the starting point
              to: the `label_column` value for the end point
              failure_flag: 0 if a path was found between the points
                and 1 if no path was found
              dist_km: The minimum aquatic distance (in kilometers)
                between the `from` and `to` locations.

        Notes
        -----
          The locations and the graph are assumed to be in a projection with
          meters as the unit of distance.
        """
        if locations.__class__.__name__ == "GeoDataFrame":
            locs = locations
        else:
            locs = gpd.read_file(locations)
        lcol = label_column

        combs = combinations(locs.index, 2)
        path_rows = []
        for indA, indB in combs:
            A = locations.loc[indA]
            B = locations.loc[indB]
            try:
                pth = self.coastal_fish_distance(A.geometry, B.geometry,
                                                 weights=weights, raise_fail=True)
                fail_flag = 0
            except NetworkXNoPath:
                pth = LineString()
                fail_flag = 1
                if raise_fail:
                    raise NetworkXNoPath("No network path between {} and {}.".format(A[lcol], B[lcol]))
            path_rows.append([A[lcol], B[lcol], fail_flag, pth])

        pathdf = gpd.GeoDataFrame(path_rows, columns=['from', 'to', 'failure_flag', 'geometry'])
        pathdf.crs = locations.crs
        pathdf['dist_km'] = pathdf.length * 1e-3
        return pathdf
