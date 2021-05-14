#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `Attribute` class."""


import unittest
from ndexncipidloader.network import Attribute
from ndex2.nice_cx_network import NiceCXNetwork


class TestAttribute(unittest.TestCase):
    """Tests for `NetworkEdgeAttribute` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor(self):
        nef = Attribute(name='name', value='value', data_type='data_type')
        self.assertEqual('name', nef.get_name())
        self.assertEqual('value', nef.get_value())
        self.assertEqual('data_type', nef.get_data_type())

    def test_str(self):
        nef = Attribute()
        self.assertEqual('name=None, values=None, type=None', str(nef))

    def test_getters_and_setters(self):
        nef = Attribute()
        self.assertEqual(None, nef.get_name())
        self.assertEqual(None, nef.get_value())
        self.assertEqual(None, nef.get_data_type())
        self.assertEqual(None, nef.get_equality_fail_reason())

        nef.set_name('name')
        nef.set_value('value')
        nef.set_data_type('data_type')
        self.assertEqual('name', nef.get_name())
        self.assertEqual('value', nef.get_value())
        self.assertEqual('data_type', nef.get_data_type())

    def test_add_attribute_to_edge(self):
        nef = Attribute(name='foo', value='someval')
        net = NiceCXNetwork()
        node_one = net.create_node('node1')
        node_two = net.create_node('node2')
        edge_id = net.create_edge(edge_source=node_one, edge_target=node_two,
                                  edge_interaction='activates')

        nef.add_attribute_to_edge(net_cx=net, edge_id=edge_id)

        e_attrs = net.get_edge_attributes(edge_id)
        self.assertEqual(1, len(e_attrs))
        self.assertEqual('foo', e_attrs[0]['n'])
        self.assertEqual('someval', e_attrs[0]['v'])

    def test_add_attribute_to_node(self):
        nef = Attribute(name='foo', value='someval')
        net = NiceCXNetwork()
        node_one = net.create_node('node1')

        nef.add_attribute_to_node(net_cx=net, node_id=node_one)

        e_attrs = net.get_node_attributes(node_one)
        self.assertEqual(1, len(e_attrs))
        self.assertEqual('foo', e_attrs[0]['n'])
        self.assertEqual('someval', e_attrs[0]['v'])

    def test__eq__(self):
        # try with self unset
        nef = Attribute()
        self.assertTrue(nef == nef)

        # try where other is none
        self.assertFalse(nef == None)
        self.assertEqual('Other is None',
                         nef.get_equality_fail_reason())

        # try where names differ
        nef1 = Attribute(name='name1', value='foo',
                                    data_type='string')
        nef2 = Attribute(name='name2', value='foo',
                                    data_type='string')
        self.assertFalse(nef1 == nef2)
        self.assertEqual('name1 name does not match name2',
                         nef1.get_equality_fail_reason())
        self.assertFalse(nef2 == nef1)
        self.assertEqual('name2 name does not match name1',
                         nef2.get_equality_fail_reason())

        # try where values differ
        nef1.set_name('name')
        nef2.set_name(nef1.get_name())
        nef1.set_value('foo')
        nef2.set_value('blah')
        self.assertFalse(nef1 == nef2)
        self.assertEqual('foo value does not match blah',
                         nef1.get_equality_fail_reason())
        self.assertFalse(nef2 == nef1)
        self.assertEqual('blah value does not match foo',
                         nef2.get_equality_fail_reason())

        # try where one value is list and other is not
        nef2.set_value(['1', '2'])
        self.assertFalse(nef1 == nef2)
        self.assertEqual('foo value does not match [\'1\', \'2\']',
                         nef1.get_equality_fail_reason())
        self.assertFalse(nef2 == nef1)
        self.assertEqual('[\'1\', \'2\'] value does not match foo',
                         nef2.get_equality_fail_reason())

        # try where both are lists but differ
        nef1.set_value(['3', '4'])
        nef2.set_value(['1', '2'])
        self.assertFalse(nef1 == nef2)
        self.assertEqual("['3', '4'] value does not match ['1', '2']",
                         nef1.get_equality_fail_reason())
        self.assertFalse(nef2 == nef1)
        self.assertEqual("['1', '2'] value does not match ['3', '4']",
                         nef2.get_equality_fail_reason())

        # try where both lists match
        nef2.set_value(nef1.get_value())

        self.assertTrue(nef1 == nef2)
        self.assertEqual(None, nef1.get_equality_fail_reason())

        self.assertTrue(nef2 == nef1)
        self.assertEqual(None, nef2.get_equality_fail_reason())

        # try where one data type is none and other is str
        nef1.set_data_type('string')
        nef2.set_data_type(None)
        self.assertTrue(nef1 == nef2)
        self.assertEqual(None, nef1.get_equality_fail_reason())

        self.assertTrue(nef2 == nef1)
        self.assertEqual(None, nef2.get_equality_fail_reason())

        # try where data types differ
        nef1.set_data_type('string')
        nef2.set_data_type('boolean')
        self.assertFalse(nef1 == nef2)
        self.assertEqual('data types differ', nef1.get_equality_fail_reason())

        self.assertFalse(nef2 == nef1)
        self.assertEqual('data types differ', nef2.get_equality_fail_reason())









