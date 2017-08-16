#!/usr/bin/env python
# -*- coding: utf-8 -*-


### Philosophy for dependencies:
###     Care more that your library successfully
###     calls the function, as opposed to the function 
###     being called every time a test is run.

"""Tests for `pyriv` package."""

import unittest
import mock
import tempfile
import os, shutil, tempfile
import doctest
from nose.tools import assert_equal, assert_not_equal, assert_raises, assert_list_equal, assert_false, assert_true, \
    assert_is_instance
from nose.tools.nontrivial import raises

import networkx as nx
from pyriv import graph_prep as GraphBuilder
from pyriv import river_graph as RiverGraph
from pyriv import point_processor as PointProcessor 


global testgraph, datapath, test_dir

class TestPyriv(unittest.TestCase):
    """Tests for `pyriv` package."""

    # NOTE: mock.patch works in the *opposite* direction
    #       arguments to the method following it!!!

    def setUp(self):
        """Set up test fixtures."""
        # make a temp directory for temp outfiles
        global test_dir, datapath, testgraph
        self.test_dir = tempfile.mkdtemp()
        self.datapath = 'testdata/'
        self.testgraph = GraphBuilder.GraphBuilder(self.datapath+'test.shp').graph
        print "setup" + str(type(self.testgraph))

    def tearDown(self):
        """Tear down test fixtures."""
        shutil.rmtree(self.test_dir)

    # def test_command_line_interface(self):
    #     """Test the CLI."""
    #     runner = CliRunner()
    #     result = runner.invoke(cli.main)
    #     assert result.exit_code == 0
    #     assert 'pyriv.cli.main' in result.output
    #     help_result = runner.invoke(cli.main, ['--help'])
    #     assert help_result.exit_code == 0
    #     assert '--help  Show this message and exit.' in help_result.output

####### --------------------                    -----------------------------    
####### -------------------- GRAPHBUILDER TESTS -----------------------------    
####### --------------------                    -----------------------------    

    # NOTE:   it may be (probably) the case that breaking the read-in process
    #         into multiple tests is redundant; consequence of being new to this

    @mock.patch('pyriv.river_graph.nx.DiGraph', autospec=True) #potentially move to setup (?)
    @mock.patch('pyriv.graph_prep.nx.read_shp', autospec=True)
    def test_create_graph_from_shp(self, mock_nx_read_shp, mock_river_graph_init):
        """Test if a networkx graph can be created from shp."""
        mock_nx_read_shp.return_value = self.testgraph 
        mock_river_graph_init.return_value = self.testgraph
      
        success = GraphBuilder.GraphBuilder(self.datapath+'test.shp').graph
        assert_list_equal(success.edges(), self.testgraph.edges())

    #   Note: this test will fail. See Issue #9
    # @mock.patch('pyriv.river_graph.nx.DiGraph', autospec=True)
    # @mock.patch('pyriv.graph_prep.nx.read_graphml', autospec=True)
    # def test_create_graph_from_graphml(self, mock_nx_read_graphml, mock_river_graph_init):
    #     """Test if a networkx graph can be created from graphml."""
    #     mock_nx_read_graphml.return_value = self.testgraph
    #     mock_river_graph_init.return_value = self.testgraph

    #     success = GraphBuilder.GraphBuilder(self.datapath+'test.graphml').graph
    #     assert_list_equal(sucess.edges(), self.testgraph.edges())

    #outgoing command; test "expect to send"
    @mock.patch('pyriv.graph_prep.nx.write_gpickle', autospec=True)
    @mock.patch('pyriv.graph_prep.nx.DiGraph', autospec=True)
    def test_graph_saving_gpickle(self, mock_nx_digraph, mock_nx_write_gpickle):
        """Test graph saving to gpickle functionality."""

        mock_nx_digraph.return_value = self.testgraph
        pass


####### --------------------                    -----------------------------    
####### --------------------   SNAPTOOL TESTS   -----------------------------    
####### --------------------                    -----------------------------
    
    #incoming command; test "(direct) public side effects"
    @mock.patch('pyriv.snapping_processor.nx.copy', autospec=True)
    @mock.path('pyriv.snapping_processor.missing_edges_list', autospec=True)
    @mock.path('pyriv.snapping_processor.add_missing_edges', autospec=True)
    def test_network_snapping_tool(self, mock_nx_copy, ):


