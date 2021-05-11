#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `NetworkNodeAttribute` class."""

import unittest
from ndexncipidloader.network import NetworkNodeAttribute
from ndex2.nice_cx_network import NiceCXNetwork


class TestNetworkNodeAttribute(unittest.TestCase):
    """Tests for `NetworkNodeAttribute` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor(self):
        nef = NetworkNodeAttribute(name='name', value='value', data_type='data_type')
        self.assertEqual('name', nef.get_name())
        self.assertEqual('value', nef.get_value())
        self.assertEqual('data_type', nef.get_data_type())

    def test_getters_and_setters(self):
        nef = NetworkNodeAttribute()
        self.assertEqual(None, nef.get_name())
        self.assertEqual(None, nef.get_value())
        self.assertEqual(None, nef.get_data_type())

        nef.set_name('name')
        nef.set_value('value')
        nef.set_data_type('data_type')
        self.assertEqual('name', nef.get_name())
        self.assertEqual('value', nef.get_value())
        self.assertEqual('data_type', nef.get_data_type())

    def test_add_attribute_to_node(self):
        nef = NetworkNodeAttribute(name='foo', value='someval')
        net = NiceCXNetwork()
        node_one = net.create_node('node1')

        nef.add_attribute_to_node(net_cx=net, node_id=node_one)

        e_attrs = net.get_node_attributes(node_one)
        self.assertEqual(1, len(e_attrs))
        self.assertEqual('foo', e_attrs[0]['n'])
        self.assertEqual('someval', e_attrs[0]['v'])
