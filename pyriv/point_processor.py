import geopandas as gpd
import pandas as pd
from river_graph import point_to_tuple
from shapely.geometry import Point

def river_distances(shp_file, riv_graph, node_distance=False):
	"""

	Parameters
	---------
	shp_file : String
		Input filepath to a .shp file that 
	riv_graph : pyriv.RiverGraph
	node_distance : {False, True}

	Returns
	-------
	pts : ??

	Notes
	-----
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

def save_points(riv_dist_gdf, outfile):
	riv_dist_gdf.drop(['path'], axis=1).to_file(outfile)

def save_paths(riv_dist_gdf, outfile):
	outgdf = riv_dist_gdf.set_geometry('path', drop=True)
	zero_ind = outdf.index[outdf.riv_dist_km==0.0]
	outdf.drop(zero_ind).to_file(outfile)