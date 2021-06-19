#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `INDRAAnnotator` class."""

import os
import logging
import unittest
from ndex2.nice_cx_network import NiceCXNetwork

from ndexncipidloader import ndexloadncipid
from ndexncipidloader.ndexloadncipid import INDRAAnnotator
from ndexncipidloader.network import NetworkEdgeFactory
from ndexncipidloader.network import NetworkEdge
from ndexncipidloader.network import Attribute


@unittest.skipIf(ndexloadncipid.INDRA_LOADED is False,
                 "ndexindraloader package not found. Skipping tests")
class TestINDRAAnnotator(unittest.TestCase):
    """Tests for `ProteinFamilyNodeMemberRemover` class."""

    def setUp(self):
        """Set up test fixtures, if any."""
        # set the logging level to debug cause some code
        # will only run in this logging level
        x = logging.getLogger('ndexncipidloader.ndexloadncipid')
        x.setLevel(logging.DEBUG)

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_description(self):
        annotator = INDRAAnnotator()
        self.assertTrue('Adds INDRA edge annotations' in annotator.get_description())

    def test_remove_edges_from_network(self):
        net = NiceCXNetwork()

        node_one = net.create_node('node1')
        node_two = net.create_node('node2')

        node_three = net.create_node('node3')
        edge_one = net.create_edge(edge_source=node_one, edge_target=node_two)
        edge_two = net.create_edge(edge_source=node_two, edge_target=node_one)

        edge_three = net.create_edge(edge_source=node_one, edge_target=node_three)

        edge_four = net.create_edge(edge_source=node_two, edge_target=node_three)

        self.assertEqual(4, len(net.get_edges()))

        net_edge_fac = NetworkEdgeFactory()
        edge_lists = [[net_edge_fac.get_network_edge_from_network(net_cx=net,
                                                                 edge_id=edge_one),
                       net_edge_fac.get_network_edge_from_network(net_cx=net,
                                                                  edge_id=edge_two)],
                      [net_edge_fac.get_network_edge_from_network(net_cx=net,
                                                                  edge_id=edge_three)]]

        annotator = INDRAAnnotator()
        annotator._remove_edges_from_network(net_cx=net, edge_lists=edge_lists)
        self.assertEqual(1, len(net.get_edges()))
        self.assertEqual(edge_four, list(net.get_edges())[0][0])




