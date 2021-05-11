#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `NetworkEdge` class."""

import unittest
from ndexncipidloader.network import NetworkEdge
from ndexncipidloader.network import NetworkEdgeAttribute
from ndexncipidloader.network import NetworkEdgeFactory
from ndex2.nice_cx_network import NiceCXNetwork


class TestNetworkEdge(unittest.TestCase):
    """Tests for `NetworkEdge` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor_and__str__(self):
        ne = NetworkEdge()
        self.assertEqual(None, ne.get_attributes())
        self.assertEqual(None, ne.get_id())
        self.assertEqual(None, ne.get_source_node_id())
        self.assertEqual(None, ne.get_target_node_id())
        self.assertEqual(None, ne.get_interaction())
        self.assertEqual(None, ne.get_source_node_name())
        self.assertEqual(None, ne.get_target_node_name())

        self.assertEqual('po=None, s=None (None), t=None (None), i=None',
                         str(ne))

        # try with values set
        ne = NetworkEdge(edge_id=1, source_node_id=2, source_node_name='src',
                         target_node_id=3, target_node_name='target',
                         interaction='interaction',
                         attributes=[])
        self.assertEqual([], ne.get_attributes())
        self.assertEqual(1, ne.get_id())
        self.assertEqual(2, ne.get_source_node_id())
        self.assertEqual(3, ne.get_target_node_id())
        self.assertEqual('interaction', ne.get_interaction())
        self.assertEqual('src', ne.get_source_node_name())
        self.assertEqual('target', ne.get_target_node_name())

        self.assertEqual('po=1, s=2 (src), t=3 (target), i=interaction',
                         str(ne))

    def test_getters_and_setters(self):
        ne = NetworkEdge()
        ne.set_attributes([])
        ne.set_id(1)
        ne.set_source_node_id(2)
        ne.set_target_node_id(3)
        ne.set_source_node_name('src')
        ne.set_target_node_name('target')
        ne.set_interaction('interaction')

        self.assertEqual([], ne.get_attributes())
        self.assertEqual(1, ne.get_id())
        self.assertEqual(2, ne.get_source_node_id())
        self.assertEqual(3, ne.get_target_node_id())
        self.assertEqual('interaction', ne.get_interaction())
        self.assertEqual('src', ne.get_source_node_name())
        self.assertEqual('target', ne.get_target_node_name())

    def test_add_edge_to_network(self):
        net = NiceCXNetwork()
        src_node = net.create_node('node1')
        target_node = net.create_node('node2')

        attrs = []

        attrs.append(NetworkEdgeAttribute(name='foo', value='ha'))
        attrs.append(NetworkEdgeAttribute(name='foo2', value='ha2',
                                          data_type='string'))

        ne = NetworkEdge(edge_id=1, interaction='activates', attributes=attrs)

        self.assertEqual(0, len(net.get_edges()))
        edge_id = ne.add_edge_to_network(net_cx=net, source_node_id=src_node,
                                         target_node_id=target_node)
        self.assertEqual(1, len(net.get_edges()))

        nef = NetworkEdgeFactory()
        new_ne = nef.get_network_edge_from_network(net_cx=net, edge_id=edge_id)
        self.assertEqual('node1', new_ne.get_source_node_name())
        self.assertEqual('node2', new_ne.get_target_node_name())
        self.assertEqual(src_node, new_ne.get_source_node_id())
        self.assertEqual(target_node, new_ne.get_target_node_id())
        self.assertEqual('activates', new_ne.get_interaction())

        self.assertEqual(2, len(new_ne.get_attributes()))
        a_dict = {}
        for attr in new_ne.get_attributes():
            a_dict[attr.get_name()] = attr

        self.assertEqual('ha', a_dict['foo'].get_value())
        self.assertEqual('ha2', a_dict['foo2'].get_value())

        new_ne.remove_edge_from_network(net_cx=net)

        self.assertEqual(0, len(net.get_edges()))

    def test_remove_edge_from_network(self):
        net = NiceCXNetwork()
        src_node = net.create_node('node1')
        target_node = net.create_node('node2')
        edge_id = net.create_edge(edge_source=src_node, edge_target=target_node)
        self.assertEqual(1, len(net.get_edges()))
        nef = NetworkEdgeFactory()
        ne = nef.get_network_edge_from_network(net_cx=net, edge_id=edge_id)

        ne.remove_edge_from_network(net_cx=net)
        self.assertEqual(0, len(net.get_edges()))
