from graph_prep import GraphBuilder
from river_graph import RiverGraph
import numpy as np
#import geopandas as gpd


def printPaths(gr, starts, g):
	paths = starts.geometry.apply(lambda g: gr.graph.shortest_path_to_coast(gr.graph.closest_node(tuple(np.array(g.coords)[0]))))
	for leng in [p.length * .001 for p in paths]:
		print "{:.2f} km".format(leng)

def main():

	gb = GraphBuilder('../../data/sasap/kusko_flowlines.shp')

	gb.graph = gb.prune_network(True)
	#gb.write_gml('kusko_prunedtry0.gml')
	gb.write_graphml('kusko_prunedtry0.graphml')
	#gb.write_gpickle('kusko_prunedtry0.pickle')

	gb_graphml = GraphBuilder('kusko_prunedtry0.graphml')
	#gb_gml = GraphBuilder('kusko_prunedtry0.gml')

	starts = gpd.read_file('../../data/sasap/test_points.shp')
	g = starts.geometry[0]
	tuple(np.array(g.coords)[0])

	printPaths(gb_graphml, starts, g)
	#printPaths(gb_gml, starts, g)



main()