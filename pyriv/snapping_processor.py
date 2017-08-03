import sys
import networkx 		as nx
import numpy 			as np
import geopandas 		as gpd
from shapely.geometry 	import Point
from multiprocessing	import Pool
from functools			import partial
from types				import *
from pyriv 				import graph_prep 	as GraphBuilder
from pyriv 				import river_graph	as river_graph



def add_node_ids_copy(G):
	"""
    Parameters
    networkx graph
    Returns copy of graph with node ID attr. added
    
    Returns a copy of the graph with node ID attr (an int) added, 
    used for finding problem nodes easier / quicker node identification
    Must use a dict with the nodes' as keys and the attributes listed
    """
    copy = G.copy()
    return add_node_ids_mutate(copy)

def node_to_dict(G):
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
    Params: networkx graph
    Return: referfence to modified graph
    
    Modifies input graph with node ID attr (an int) added, 
    used for finding problem nodes easier / quicker node identification
    Must use a dict with the nodes' as keys and the attributes listed
    """
    attr_dict = {}
    i = 0
    keys = G.nodes()
    for node in keys:
        attr_dict[node] = i
        i = i+1
    nx.set_node_attributes(G, 'ID', attr_dict)
    return G

def deadends_old(self):
    """
    For a directed graph, find the nodes that have out_degree == 0
    Note: this is the older, slower method
    """
    thedict = {}
    # get a dictionary of { node: out_degree }
    degdict = self.out_degree(self.nodes())
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


def deadend_coords_to_keys(G, ends):
    assert type(ends) is TupleType, "collection of end nodes is not a tuple"
    ends = list(ends)
    
    
    alldict = nx.get_node_attributes(G, 'ID')
    
    keylist = []
    for i in range(0, len(ends)):
        keylist.append(alldict[ends[i]])
        
    return keylist



def find_missing_edges_par(ends, nodes, threshold, numproc):
    """
    Params: ends = deadends found dict, allnodes = iterable of all nodes; threshold = dist. thresh
    Returns: a list of missing edges to add where dist < thresh. b/w the deadend & other node
    
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
    inner function for || version of finding missing edges
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
    Params: G = graph to add edges to; list of edges to add w/ each entry an [x,y] list
    Returns: nothing
    
    IMPORTANT: MUTATES the graph passed in
    """
    for coord in edges_to_add:
        G.add_edge(allnodes[coord[0]], allnodes[coord[1]])


def missing_edges_to_shp(pts, filepath, projection_code):
    """
    Params: nodes/pts describing missing edges; relative filepath .shp to write to; projection_code to convert
            coords (should match the maps you started with for overlaying purposes)
    Returns: none (writes to output file) 
    
    
    """
    pts = list(pts.values())
    kcdf = gpd.GeoDataFrame({'geometry': [Point(n) for n in pts]})
    kcdf.crs = {'init': projection_code}
    kcdf.to_file(filepath)


def missing_edges_list(graph_to_snap, dist_thresh, outfile):
	copy_w_id = add_node_ids_copy(graph_to_snap.graph)
	ends = copy_w_id.deadends()
	ends_list = deadend_coords_to_keys(ends)
	allnodesdict = node_to_dict(copy_w_id)

	# default number of processes set to 8 
	num_proc = 8
	missing_edges = find_missing_edges_par(ends_list, allnodesdict, dist_thresh, num_proc)
	return allnodesdict, missing_edges

def snap_graph(graph_to_snap, dist_thresh, outfile):
	graph_copy = graph_to_snap.copy()
	allnodesdict, missing_edges = missing_edges_list(graph_to_snap, dist_thresh, outfile)
	add_missing_edges(graph_copy, missing_edges, allnodesdict)
	return graph_copy

############## ---- STAT METHODS -----------------


def component_stats(G, verbose):
    """
    Params: G = graph to collect stats on, verbose flag
    Returns: nothing
    
    Prints out various relevent stats about graphs concerning components
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
    """
    Print general stats about graph G
    Can add more as/if needed...
    """
    print "Number of edges: "+ str(G.number_of_edges())
    print "Number of nodes: "+str(G.number_of_nodes()) 

############## ---- UNUSED / LEGACY METHODS -----------------


def find_missing_edges_seq(G, ends, allnodes, threshold):
    """
    Params:
    Returns: list of missing edges

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