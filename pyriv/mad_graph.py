import networkx as nx
import geopandas as gpd
from shapely.geometry import LineString

class MadGraph(nx.Graph):
    """
    This class represents a network graph of routes that can be traversed without
    crossing land to find Minimum Aquatic Distance (MAD).
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
