#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `NetworkNode` class."""

import unittest
from ndexncipidloader.network import NetworkNode
from ndexncipidloader.network import NetworkNodeFactory
from ndex2.nice_cx_network import NiceCXNetwork


class TestNetworkNode(unittest.TestCase):
    """Tests for `NetworkNode` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor_and__str__(self):
        nn = NetworkNode()
        self.assertEqual(None, nn.get_attributes())
        self.assertEqual(None, nn.get_name())
        self.assertEqual(None, nn.get_id())
        self.assertEqual(None, nn.get_represents())
        self.assertEqual('@id=None, name=None, represents=None',
                         str(nn))

        # try with values set
        nn = NetworkNode(node_id=1, name='name',
                         represents='represents', attributes=[],
                         network_edge_factory='ha')
        self.assertEqual([], nn.get_attributes())
        self.assertEqual('name', nn.get_name())
        self.assertEqual(1, nn.get_id())
        self.assertEqual('represents', nn.get_represents())
        self.assertEqual('@id=1, name=name, represents=represents',
                         str(nn))
        # digging into class here, but easier to test
        self.assertEqual('ha', nn._nef)

    def test_getters_and_setters(self):
        nn = NetworkNode(node_id=1, name='name',
                         represents='represents', attributes=[])

        nn.set_name('name2')
        nn.set_represents('represents2')
        nn.set_id(2)
        nn.set_attributes(['ha'])
        self.assertEqual('name2', nn.get_name())
        self.assertEqual('represents2', nn.get_represents())
        self.assertEqual(2, nn.get_id())
        self.assertEqual(['ha'], nn.get_attributes())

    def test_remove_node_from_network_no_edges_no_attrs(self):
        net = NiceCXNetwork()
        node_id = net.create_node('node1', node_represents='ha')
        nef = NetworkNodeFactory()
        nn = nef.get_network_node_from_network(net_cx=net, node_id=node_id)

        nn.remove_node_from_network(net_cx=net)

        self.assertEqual(0, len(net.get_nodes()))

    def test_remove_node_from_network_with_edges_and_attrs(self):
        net = NiceCXNetwork()
        node_id = net.create_node('node1', node_represents='ha')
        other_node_id = net.create_node('node2', node_represents='ha')
        edge_id = net.create_edge(edge_source=node_id, edge_target=other_node_id,
                                  edge_interaction='activates')

        net.set_node_attribute(node_id, 'nodeattr', values='foo')
        net.set_node_attribute(node_id, 'nodeattr2', values='foo2')

        net.set_edge_attribute(edge_id, 'edgeattr', values='hello')
        net.set_edge_attribute(edge_id, 'edgeattr2', values='hello2')

        nef = NetworkNodeFactory()
        nn = nef.get_network_node_from_network(net_cx=net, node_id=node_id)

        nn.remove_node_from_network(net_cx=net)

        self.assertEqual(1, len(net.get_nodes()))
        self.assertEqual(0, len(net.get_edges()))
