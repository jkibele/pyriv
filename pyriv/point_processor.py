import geopandas as gpd
import pandas as pd
from river_graph import point_to_tuple
from shapely.geometry import Point

def river_distances(shp_file, riv_graph, node_distance=False):
	"""
	Assumes projection units are meters for the moment.
	"""
	pnts = gpd.read_file(shp_file)
	# function to find closest node for a point geometry
	cl_nd = lambda p: riv_graph.closest_node(point_to_tuple(p))
	path_find = lambda p: riv_graph.shortest_path_to_coast(cl_nd(p))
	pnts['path'] = pnts.geometry.apply(path_find)
	pth_dist = lambda p: p.length * 0.001
	pnts['riv_dist_km'] = pnts.path.apply(pth_dist)
	if node_distance:
		pnts['nearest_node'] = pnts.geometry.apply(lambda p: Point(cl_nd(p)))
		d = {}
		for i, row in pnts.iterrows():
			d[i] = row['geometry'].distance(row['nearest_node'])
		dser = pd.Series(d, name='node_dist')
		pnts = pnts.join(dser)
	return pnts
