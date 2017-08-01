import os
import numpy as np
import networkx as nx
from river_graph import RiverGraph

class GraphBuilder(object):

    def __init__(self, file_path, coastline_shp=None, calc_dist_weights=True):
        """
        Build a graph from a shapefile and provide methods to prune, or read
        a gpickle file and convert it to a `RiverGraph`.
        """
        if os.path.exists(file_path) == True and os.path.isfile(file_path) == True:
            self.coast_fn = coastline_shp
            if os.path.splitext(file_path)[1] == '.shp':
                g = nx.read_shp(file_path)
                print(g)
                print type(g)
                print(g.edges())
            elif os.path.splitext(file_path)[1] == '.graphml':
                g = nx.read_graphml(file_path)
            elif os.path.splitext(file_path)[1] == '.gpickle':
                g = nx.read_gpickle(file_path)
            else:
                #assert False, "path does not have a valid extention (.shp, .graphml, .gpickle): %s" % file_path
                print "path does not have a valid extension (.shp, .graphml, .gpickle)"

            self.graph = RiverGraph(data=g, coastline_shp=self.coast_fn)
            if calc_dist_weights:
                print "Weighting Edges with Distances"
                self.graph = self.graph.weight_edges()
            #return True
        else:
            if os.path.exists(file_path) == False:
                print "file does not exist"
            if os.path.isfile(file_path) == False:
                print "path is not a file or does not exist: %s" % file_path
            #return False
        

    def prune_network(self, verbose=False):
        """
        Remove subgraphs of the network that do not connect to the coastline.
        """
        # weakly_connected_component_subgraphs doesn't work with the 
        # RiverGraph subclass, so we have to switch it back to DiGraph
        sg = self.graph
        coast_n = 0
        noncoast_n = 0
        g_list = []
        cps = nx.weakly_connected_component_subgraphs(sg)
        for cp in cps:
            if has_rivermouth(cp.nodes(), sg):
                g_list.append(cp)
                coast_n += 1
            else:
                noncoast_n += 1
        if verbose:
            print "{} graphs with coastal nodes, {} without.".format(coast_n, noncoast_n)
        return RiverGraph(data=nx.compose_all(g_list), coastline_shp=self.coast_fn)

    def write_gpickle(self, out_file_path):
        """
        This writes the DiGraph to a gpickle.
        """
        nxg = nx.DiGraph(data=self.graph)
        nx.write_gpickle(nxg, out_file_path)

    def write_graphml(self, out_file_path):
        """
        This writes the DiGraph to a GraphML file.
        """
        nxg = nx.DiGraph(data=self.graph)
        nx.write_graphml(nxg, out_file_path)


def has_rivermouth(node_list, sg):
    """
    Given a list of nodes, return `True` if at least one node is
    coastal. Otherwise, return `False`. This is essentially the same as
    `RiverGraph.has_coast_node`, but including this function here stops me
    from having to convert subgraphs to RiverGraphs when I prune the network.
    """
    return np.apply_along_axis(sg.is_rivermouth, 1, np.array(node_list)).any()
