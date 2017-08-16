import errno
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
        else:
            if os.path.exists(file_path) == False:
                print "file does not exist"
            if os.path.isfile(file_path) == False:
                print "path is not a file or does not exist: %s" % file_path
        

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


    def splitpath(self, filename):
        """
        Splits full path into path and filename, respectively.
        """
        base = os.path.basename(filename)
        tmp = filename.split(base)
        return (tmp[0], base)

    def validate_file(self, filepath, mode='read'):
        """
        This method checks the validity of a file path & name
        """
        retval=()
        stop = True
        while(stop):
            patharr = self.splitpath(filepath)
            if mode == 'read':
                if os.path.exists(filepath) == True and os.path.isfile(filepath) == True:
                    retval = retval + (patharr[0], patharr[1])
                    stop = False
                else:
                    filepath = raw_input("You entered and invalid filename. Please try again:") #either not a file or does not exist in that directory
            elif mode == 'write':
                if os.path.exists(filepath) == False:
                    retval = retval + (patharr[0], patharr[1])
                    stop = False
                else:
                    filepath = raw_input("That filename already exists. Please try again:")
        return retval

    def write_shp(self, out_file_path):
        """
        This writes the DiGraph to GIS .shp project files.

        Note 1: You will need to set the CRS manually in a GIS program;
                that is not in this code as the processing can be too costly. 
        Note 2: Will only create one level of non-existent directories, else an OSError.

        """
        #arr = self.splitpath(out_file_path)
        #tmpstr = arr[0] + os.path.splitext(arr[1])[0]
        #os.mkdir(tmpstr)
        (path, filename) = self.validate_file(out_file_path, mode='write')
        splitname = os.path.splitext(filename)

        #create a directory with that path to avoid overwriting edges.shp
        #when writing multiple networks out to the same directory
        try:
            os.mkdir(path + splitname[0])
            print "made directory"
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        fullpath = path + splitname[0] + '/' + filename
        nx.write_shp(self.graph, fullpath)

def has_rivermouth(node_list, sg):
    """
    Given a list of nodes, return `True` if at least one node is
    coastal. Otherwise, return `False`. This is essentially the same as
    `RiverGraph.has_coast_node`, but including this function here stops me
    from having to convert subgraphs to RiverGraphs when I prune the network.
    """
    return np.apply_along_axis(sg.is_rivermouth, 1, np.array(node_list)).any()
