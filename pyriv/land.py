import geopandas as gpd
from .common import *

from shapely.geometry import Polygon

class CoastLineGraph(nx.Graph):
	"""
	This class can turn a LineString representation of a coastline into a 
	polygon representation.
	"""
    def __init__(self, *args, **kwargs):
        """
        Build a CoastLine object.

        Parameters
        ----------

        Returns
        -------
          A CoastLine object
        """
        self = super(CoastLine, self).__init__(*args, **kwargs)
        
    @classmethod
    def read_shp(cls, shp_fn):
        dig = nx.read_shp(shp_fn, simplify=False)
        return cls(dig)
    
    def connected_subgraphs(self):
        return nx.connected_component_subgraphs(self)
    
    def rings(self):
        rings = [list(nx.dfs_preorder_nodes(sg)) for sg in self.connected_subgraphs()]
        return rings
    
    def polygons(self):
        return [Polygon(r) for r in self.rings()]
    
    def poly_geodataframe(self):
        return gpd.GeoDataFrame({'geometry': self.polygons()})

class Land(gpd.GeoDataFrame):
	"""
	This class will represent land geometry in order to generate a coastal graph.
	"""
	def __init__(self, *args, **kwargs):
	    """
	    Build a RiverGraph object.
	    
	    Parameters
	    ----------
	      
	        
	    Returns
	    -------
	      A Land object
	    """
	    self = super(Land, self).__init__(*args, **kwargs)