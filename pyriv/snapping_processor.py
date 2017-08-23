import networkx as nx
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from multiprocessing import Pool
from functools import partial
from types import *
from river_graph import RiverGraph
from graph_prep import GraphBuilder



def add_node_ids_copy(G):
	"""
	Copies graph and adds numerical ID to nodes as attribute to the copy.

	Parameters
	----------
	  G : networkx digraph

	Returns
	-------
	  mutated_copy : networkx digraph Returns copy of graph with node ID attr. added
	    copy of G with ID attribute (an int) on each node
	
	Notes: Used for finding problem nodes easier / quicker node identification
	??* Must use a dict with the nodes' as keys and the attributes listed
	"""

	assert isinstance(G, nx.DiGraph) is True, "graph is not a digraph"

	#print type(G)
	mutated_copy = add_node_ids_mutate(G.copy())
	return mutated_copy

def node_to_dict(G):
	"""
	Puts graph nodes into a dict with numerical IDs as keys. 

	Parameters
	----------
	  G : networkx digraph

	Returns
	-------
	  temp : dict, int --> node
	    A dict with int IDs mapped to nodes from graph G

	Note: Creates a way to reference nodes in a predictable, easily-indexible format
	"""
	assert isinstance(G, nx.DiGraph) is True, "graph is not a digraph"

	nodes = nx.nodes(G)
	#check if they have id first?
	#if nodes[0]
	temp = {}
	i = 0
	for node in nodes:
		temp[i] = node
		i = i+1
		
	return temp


def add_node_ids_mutate(G):
	"""
	Mutates graph and adds numerical ID to nodes as attribute to the copy.

	Parameters
	----------
	  G : networkx digraph

	Returns
	-------
	  G : networkx digraph
	    Refernce to modified graph with numerical IDs added as node attributes

	??Must use a dict with the nodes' as keys and the attributes listed
	"""
	assert isinstance(G, nx.DiGraph) is True, "graph is not a digraph"
	assert isinstance(G, dict) is False, "this should be a graph not a dict"
	attr_dict = {}
	i = 0
	keys = G.nodes()
	for node in keys:
		attr_dict[node] = i
		i = i+1
	nx.set_node_attributes(G, 'ID', attr_dict)
	return G


def deadend_coords_to_keys(G, ends):
	"""
	Parameters
	----------
	  G : networkx graph
	    Must have node attribute 'ID' of type int
	  ends : tuple of ??

	Returns
	-------
	  keylist : list 
	"""
	assert type(ends) is TupleType, "collection of end nodes is not a tuple"
	ends = list(ends)
	
	
	alldict = nx.get_node_attributes(G, 'ID')
	
	keylist = []
	for i in range(0, len(ends)):
		try:
			keylist.append(alldict[ends[i]])
		except KeyError:
			print alldict[ends[i]]

		
	return keylist



def find_missing_edges_par(ends, nodes, threshold, numproc):
	"""
	Parallel wrapper method for finding missing edges algorithm.

	Parameters
	----------
	  ends : dict of ??
	    All the deadends found
	  allnodes : iterable of all nodes
	  threshold : float
	    Max distance to snap
	  numproc: int
	    The number of processes you want to run in order to complete the task

	Returns
	-------
	  cleaned_results : list of edges
	    Missing edges to add where dist < threshold bewtween the deadend and other node 

	
	||'zd method to find missing edges (a CPU-bound task w/ no communication needs)
	
	Note 1: in order to pass multiple args to multiproc.'g map fcn. we use a partial 
																( req. Python >=2.7)
	Note 2: Pool.map maps a fcn over sequence (str, unicode, list, tuple, buffer, or xrange) so
						need to convert dict to that first
	"""
	assert type(ends) is ListType, "collection of end nodes is not a list"
	assert type(ends[0]) is IntType, "nodes in ends list arent referred to by their ID"
	
	assert type(nodes) is DictType, "collection of allnodes is not a dictionary"
	
	#need to convert the dict of dead ends to the correct format, which is:
	#    key: (key, (coord_x, coord_y))

	realdict = {}
	for key, val in nodes.iteritems():
		realdict[key]=(key, val)
		
	#assert type(allnodes[0]) is IntType, "first column of allnodes dictionary is not a node ID"
	#assert type(allnodes[0][0]) is TupleType, "second column of allnodes dictionary is not a coord tuple"
	
	
	if numproc < 2:
		print("Please request to run more than one concurrent process to avoid redundant work. If you wish to use "+
				+"the sequential version of this method, please do so directly.")
	else:
		
		
		  
		pool = Pool(numproc)
	
		result_list = pool.map(partial(find_missing_edges_par_inner, allnodes=realdict, th=threshold), ends)
		pool.close()
		pool.join()
		
		cleaned_results = []
		#clean result list before returning or else it will have empty entries
		for entry in result_list:
			if not not entry:
				if len(entry) > 1:
					for innerentry in entry:
						cleaned_results.append(innerentry)
				else:
					cleaned_results.append(entry[0])
				
		return cleaned_results 


def find_missing_edges_par_inner(key1, allnodes, th):
	"""
	Parallel inner method/algorithm for finding missing edges.

	Parameters
	----------
	  key1 : ??
	  allnodes : iterable of nodes
	    ??
	  th : float
	    Max distance to snap

	Returns
	-------
	  edges_to_add: ??
	    ??
	"""
	
	edges_to_add = []
	node1 = allnodes[key1][1]
	for i in range(0, len(allnodes)):
		key2 = allnodes[i][0]
		if key1 != key2:
			node2 = allnodes[i][1]
			xd = node1[0] - node2[0]
			yd = node1[1] - node2[1]
			dist = np.sqrt((xd*xd + yd*yd))
			if key1 != key2 and dist < th and dist > 0:
				edges_to_add.append([key1,key2])
				break
				
	return edges_to_add


def add_missing_edges(G, edges_to_add, allnodes):
	"""
	Add missing edges to a graph.

	Parameters
	----------
	  G : networkx graph
	    Graph to add edges ti
	  edges_to_add : ?? list of edges to add w/ each entry an [x,y] list
	  allnodes : ??

	Returns
	-------
	
	Note: MUTATES the graph passed in
	"""
	for coord in edges_to_add:
		G.add_edge(allnodes[coord[0]], allnodes[coord[1]])


def missing_edges_to_shp(pts, filepath, projection_code):
	"""
	Parameters
	----------
	  pts: ??
	    nodes/pts describing missing edges
	  filepath : string
	    Relative filepath .shp to write to
	  projection_code: string
	    Projection_code to convert coords (should match the maps you started with for overlaying purposes)

	Returns
	-------
	
	Note: Writes out a file.
	"""
	pts = list(pts.values())
	kcdf = gpd.GeoDataFrame({'geometry': [Point(n) for n in pts]})
	kcdf.crs = {'init': projection_code} #AA == 'epsg:3338'
	kcdf.to_file(filepath)


def missing_edges_list(graph_to_snap, dist_thresh, outfile):
	"""
	Parameters
	----------
	  * see snap_graph method

	Returns
	-------
	"""

	assert isinstance(graph_to_snap, nx.DiGraph) is True, "graph is not a digraph"
	
	copy_w_id = add_node_ids_copy(graph_to_snap)
	ends = copy_w_id.deadends()
	ends_list = deadend_coords_to_keys(copy_w_id, ends)
	allnodesdict = node_to_dict(copy_w_id)

	# default number of processes set to 8 
	num_proc = 8
	missing_edges = find_missing_edges_par(ends_list, allnodesdict, dist_thresh, num_proc)
	return (allnodesdict, missing_edges)


def snapped_graph(graph_to_snap, dist_thresh, outfile):
	"""
	Returns a snapped copy of a graph.

	Parameters
	----------
	  graph_to_snap : networkx DiGraph
	  dist_thresh : float
	    snap gaps in network <= dist_thresh
	  outfile : string
	    file to write new network to (unused currently)

	Returns
	-------
	  graph_copy : networkx DiGraph
	    copy of graph_to_snap, with gaps <= dist_thresh snapped
	"""

	try:
		graph_copy = graph_to_snap.copy()
		allnodes_missingedges = missing_edges_list(graph_copy, dist_thresh, outfile)
		add_missing_edges(graph_copy, allnodes_missingedges[1], allnodes_missingedges[0])
		return graph_copy
	except TypeError, AttributeError:
		if not isinstance(dist_thresh, float):
			print "dist_thresh is not a float/cannot be widened into a float"
		if not isinstance(outfile, str):
			print "outfile name is not a string"


############## ---- STAT METHODS -----------------


def component_stats(G, verbose):
	"""Prints out various relevent stats about graphs concerning components.

	Parameters
	----------
	  G : networkx DiGraph
	  verbose : bool
	    Set to True if you want explanations of stats

	Returns
	-------

	Note: Writes to terminal.
	"""

	explans = {}
	if verbose == True:
		explans['weakly-connected'] = "(There is an undirected path between each pair of nodes in the directed graph)"
		explans['strongly-connected'] = "(There is a directed path between each pair of nodes in the directed graph)"
		explans['semiconnected'] = ""
	else:
		explans['weakly-connected'] = ""
		explans['strongly-connected'] = ""
		explans['semiconnected'] = ""
		
	
	print "Is the graph weakly connected "+explans['weakly-connected'] +"? "+ str(nx.is_weakly_connected(G)) 
	print "Number of weakly connected components: " + str(nx.number_weakly_connected_components(G))
	print "Is the graph semiconnected "+explans['semiconnected']+ "? " + str(nx.is_semiconnected(G))
	print "Is the graph strongly connected "+explans['strongly-connected']+ "? "+ str(nx.is_strongly_connected(G))
	
	
def general_stats(G):
	"""Print general stats about graph G

	Parameters
	----------
	  G : networkx DiGraph

	Returns
	-------

	Notes
	-----
	Can add more as/if needed...
	"""

	print "Number of edges: "+ str(G.number_of_edges())
	print "Number of nodes: "+str(G.number_of_nodes()) 

############## ---- UNUSED / LEGACY METHODS -----------------


def find_missing_edges_seq(G, ends, allnodes, threshold):
	"""
	Parameters
	----------

	Returns
	-------
	 list of missing edges

	Sequential method to find missing edges. VERY VERY SLOW; use ||'zd version
	"""
	for key1, node1 in deadends.iteritems():
		for key2, node2 in allnodes.iteritems():
			#print node1
			xd = node1[0] - node2[0]
			yd = node1[1] - node2[1]
			dist = np.sqrt((xd*xd + yd*yd))
			
			if key1 != key2 and dist < threshold and dist > 0:
				graph.add_edge(key1, key2)
				print("added edge from "+str(key1)+" : "+str(node1)+" to "+ str(key2)+" : "+str(node2))
				problem_points[str(cnt)+"_A"] =key1
				problem_points[str(cnt)+"_B"] =key2
				break

def deadends_old(G):
	"""
	Parameters
	----------

	Returns
	-------

	For a directed graph, find the nodes that have out_degree == 0
	Note: this is the older, slower method
	"""

	assert isinstance(G, nx.DiGraph) is True, "graph is not a digraph"

	thedict = {}
	# get a dictionary of { node: out_degree }
	degdict = G.out_degree(G.nodes())
	# convert to a 2 column array of node and out_degree
	degarr = np.array(degdict.items())
	# get a boolean index of rows where out_degree==0
	dead_row_ind = (degarr[:,1]==0)

	
	for i in range(0, len(degarr)):
		if dead_row_ind[i] == True:
			thedict[allnodes[i][0]] = tuple(degarr[:,0][i])
			
	# use that index to get the nodes
	#dead_nodes = degarr[:,0][dead_row_ind]
	return thedict

############## ---- CLASS DEFINITION -----------------
############## ---- CLASS DEFINITION -----------------

class SnapTool(object):
    """
    A graph representation of a river network.
    """
    def __init__(self):
    	pass
