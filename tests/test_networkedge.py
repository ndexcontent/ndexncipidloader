#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `NetworkEdge` class."""

import unittest
from ndexncipidloader.network import NetworkEdge
from ndexncipidloader.network import Attribute
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

        attrs.append(Attribute(name='foo', value='ha'))
        attrs.append(Attribute(name='foo2', value='ha2',
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

    def test_add_edge_to_network_using_src_target_node_ids(self):
        net = NiceCXNetwork()
        src_node = net.create_node('node1')
        target_node = net.create_node('node2')

        ne = NetworkEdge(edge_id=1, source_node_id=src_node,
                         target_node_id=target_node,
                         interaction='activates')

        self.assertEqual(0, len(net.get_edges()))
        edge_id = ne.add_edge_to_network(net_cx=net)
        self.assertEqual(1, len(net.get_edges()))

        nef = NetworkEdgeFactory()
        new_ne = nef.get_network_edge_from_network(net_cx=net, edge_id=edge_id)
        self.assertEqual('node1', new_ne.get_source_node_name())
        self.assertEqual('node2', new_ne.get_target_node_name())
        self.assertEqual(src_node, new_ne.get_source_node_id())
        self.assertEqual(target_node, new_ne.get_target_node_id())
        self.assertEqual('activates', new_ne.get_interaction())

        self.assertEqual([], new_ne.get_attributes())
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

    def test_remove_identical_attributes_no_attrs(self):
        ne = NetworkEdge()
        # test where there are no attributes
        self.assertEqual([], ne.remove_identical_attributes())

    def test_remove_identical_attributes_one_attr(self):
        ne = NetworkEdge()
        # test where there is 1 attribute
        attrs = [Attribute(name='foo', value='bar')]
        ne.set_attributes(attrs)
        self.assertEqual([], ne.remove_identical_attributes())
        self.assertEqual(1, len(ne.get_attributes()))

    def test_remove_identical_attributes_two_diff_attrs(self):
        ne = NetworkEdge()
        # test where there is 1 attribute
        attrs = [Attribute(name='foo', value='bar'),
                 Attribute(name='blah', value='yo')]
        ne.set_attributes(attrs)
        self.assertEqual([], ne.remove_identical_attributes())
        self.assertEqual(2, len(ne.get_attributes()))

    def test_remove_identical_attributes_four_attrs_two_same(self):
        ne = NetworkEdge()
        # test where there is 1 attribute
        attrs = [Attribute(name='foo', value='bar'),
                 Attribute(name='blah', value='yo'),
                 Attribute(name='blah', value='yo'),
                 Attribute(name='ha', value='well')]
        ne.set_attributes(attrs)
        merged = ne.remove_identical_attributes()
        self.assertEqual(1, len(merged))
        self.assertEqual('blah', merged[0].get_name())
        self.assertEqual(3, len(ne.get_attributes()))

    def test_remove_identical_attributes_four_attrs_three_same(self):
        ne = NetworkEdge()
        # test where there is 4 attributes
        attrs = [Attribute(name='foo', value='bar'),
                 Attribute(name='blah', value='yo'),
                 Attribute(name='blah', value='yo'),
                 Attribute(name='blah', value='yo')]
        ne.set_attributes(attrs)
        merged = ne.remove_identical_attributes()
        self.assertEqual(2, len(merged))
        self.assertEqual('blah', merged[0].get_name())
        self.assertEqual('blah', merged[1].get_name())
        self.assertEqual(2, len(ne.get_attributes()))

    def test_remove_identical_attributes_five_attrs_two_pair_same(self):
        ne = NetworkEdge()
        # test where there is 5 attribute
        attrs = [Attribute(name='x', value='xval'),
                 Attribute(name='a', value='bar'),
                 Attribute(name='blah', value='yo'),
                 Attribute(name='a', value='bar'),
                 Attribute(name='blah', value='yo')]
        ne.set_attributes(attrs)
        merged = ne.remove_identical_attributes()
        self.assertEqual(2, len(merged))
        merged_anames = set(x.get_name() for x in merged)
        self.assertTrue('a' in merged_anames)
        self.assertTrue('blah' in merged_anames)
        self.assertEqual(2, len(merged_anames))
        self.assertEqual(3, len(ne.get_attributes()))

    def test_get_data_type_for_attributes_single_type(self):
        ne = NetworkEdge()

        # single None
        attrs = [Attribute(data_type=None)]
        self.assertEqual('list_of_string',
                         ne._get_data_type_for_attributes(attrs))

        # multiple None
        attrs = [Attribute(data_type=None),
                 Attribute(data_type=None)]
        self.assertEqual('list_of_string',
                         ne._get_data_type_for_attributes(attrs))

        # single other types
        for a in ['string', 'boolean',
                  'double', 'integer', 'long']:
            attrs = [Attribute(data_type=a)]
            self.assertEqual('list_of_' + a,
                             ne._get_data_type_for_attributes(attrs))

        # multiple of matching types
        for a in ['string', 'boolean',
                  'double', 'integer', 'long',
                  'list_of_boolean', 'list_of_double',
                  'list_of_integer', 'list_of_long',
                  'list_of_string']:
            attrs = [Attribute(data_type=a),
                     Attribute(data_type=a),
                     Attribute(data_type=a)]
            if a.startswith('list_of'):
                checkval = a
            else:
                checkval = 'list_of_' + a
            self.assertEqual(checkval,
                             ne._get_data_type_for_attributes(attrs))

    def test_get_data_type_for_attributes_two_types(self):
        ne = NetworkEdge()

        # try various combos of <type> and None
        for a in [(None, 'string'),
                  ('string', None),
                  (None, 'list_of_string'),
                  ('list_of_string', None),
                  ('boolean', None),
                  (None, 'boolean'),
                  ('list_of_integer', None),
                  (None, 'list_of_double')]:

            attrs = [Attribute(data_type=a[0]),
                     Attribute(data_type=a[1])]
            self.assertEqual('list_of_string',
                             ne._get_data_type_for_attributes(attrs))

        # try the other types
        for a in [('string', 'list_of_string'),
                  ('list_of_string', 'string'),
                  ('boolean', 'list_of_boolean'),
                  ('list_of_boolean', 'boolean'),
                  ('double', 'list_of_double'),
                  ('list_of_double', 'double'),
                  ('integer', 'list_of_integer'),
                  ('list_of_integer', 'integer'),
                  ('long', 'list_of_long'),
                  ('list_of_long', 'long')]:
            attrs = [Attribute(data_type=a[0]),
                     Attribute(data_type=a[1])]
            if a[0].startswith('list_of'):
                checkval = a[0]
            else:
                checkval = 'list_of_' + a[0]
            self.assertEqual(checkval,
                             ne._get_data_type_for_attributes(attrs))

    def test_merge_attributes_test_where_no_dtype(self):
        ne = NetworkEdge()

        attrs = [Attribute(name='x', value='1'),
                 Attribute(name='x', value='2')]

        res = ne._merge_attributes(attrs_list=attrs)
        self.assertEqual('x', res.get_name())
        self.assertEqual(2, len(res.get_value()))
        self.assertTrue('1' in res.get_value())
        self.assertTrue('2' in res.get_value())
        self.assertEqual('list_of_string', res.get_data_type())

    def test_merge_attributes_test_where_dtype_mixed_string_none(self):
        ne = NetworkEdge()

        attrs = [Attribute(name='x', value='1', data_type='string'),
                 Attribute(name='x', value='2')]

        res = ne._merge_attributes(attrs_list=attrs)
        self.assertEqual('x', res.get_name())
        self.assertEqual(2, len(res.get_value()))
        self.assertTrue('1' in res.get_value())
        self.assertTrue('2' in res.get_value())
        self.assertEqual('list_of_string', res.get_data_type())

    def test_merge_attributes_test_with_integer_dtype(self):
        ne = NetworkEdge()

        attrs = [Attribute(name='x', value=1, data_type='integer'),
                 Attribute(name='x', value=[2, 3],
                           data_type='list_of_integer')]

        res = ne._merge_attributes(attrs_list=attrs)
        self.assertEqual('x', res.get_name())
        self.assertEqual(3, len(res.get_value()))
        self.assertTrue(1 in res.get_value())
        self.assertTrue(2 in res.get_value())
        self.assertTrue(3 in res.get_value())
        self.assertEqual('list_of_integer', res.get_data_type())

    def test_merge_attributes_test_with_str_dtype(self):
        ne = NetworkEdge()

        attrs = [Attribute(name='x', value=['1'],
                           data_type='list_of_string'),
                 Attribute(name='x', value='3')]

        res = ne._merge_attributes(attrs_list=attrs)
        self.assertEqual('x', res.get_name())
        self.assertEqual(2, len(res.get_value()))
        self.assertTrue('1' in res.get_value())
        self.assertTrue('3' in res.get_value())
        self.assertEqual('list_of_string', res.get_data_type())

    def test_merge_duplicate_attributes(self):

        # test merge with None for attributes
        ne = NetworkEdge()
        ne.merge_duplicate_attributes()

        # test merge with only 1 attribute
        ne.set_attributes([Attribute(name='foo', value='bar')])
        ne.merge_duplicate_attributes()
        self.assertEqual(1, len(ne.get_attributes()))
        self.assertEqual('foo', ne.get_attributes()[0].get_name())

        # test merge with two identical attributes
        ne.set_attributes([Attribute(name='foo', value='bar'),
                           Attribute(name='foo', value='bar')])

        ne.merge_duplicate_attributes()
        self.assertEqual(1, len(ne.get_attributes()))
        self.assertEqual('foo', ne.get_attributes()[0].get_name())

        # test merge with two identical and 1 different
        ne.set_attributes([Attribute(name='foo', value='bar'),
                           Attribute(name='foo', value='bar'),
                           Attribute(name='other', value='yo')])

        ne.merge_duplicate_attributes()
        self.assertEqual(2, len(ne.get_attributes()))
        a_names = set(x.get_name() for x in ne.get_attributes())
        self.assertTrue('foo' in a_names)
        self.assertTrue('other' in a_names)

    def test_merge(self):

        # test against None
        ne = NetworkEdge()
        ne.merge_edge(None)

        # test common scenario
        ne_one = NetworkEdge(attributes=[Attribute(name='directed', value=True,
                                                   data_type='boolean')])
        ne_two = NetworkEdge(attributes=[Attribute(name='directed', value=True,
                                                   data_type='boolean')])

        ne.merge_edge(ne_two)
        self.assertTrue(1, len(ne_one.get_attributes()))
        self.assertEqual('directed', ne.get_attributes()[0].get_name())
        self.assertEqual(True, ne.get_attributes()[0].get_value())
        self.assertEqual('boolean', ne.get_attributes()[0].get_data_type())

        # test another common scenario
        ne_one = NetworkEdge(attributes=[Attribute(name='citation',
                                                   value=['1', '2'],
                                                   data_type='list_of_string'),
                                         Attribute(name='directed', value=True,
                                                   data_type='boolean')
                                         ])
        ne_two = NetworkEdge(attributes=[Attribute(name='citation',
                                                   value=['3', '1', '2'],
                                                   data_type='list_of_string'),
                                         Attribute(name='directed', value=True,
                                                   data_type='boolean')
                                         ])

        ne.merge_edge(ne_two)
        self.assertEqual(2, len(ne.get_attributes()))
        a_names = set(x.get_name() for x in ne.get_attributes())
        self.assertTrue('citation' in a_names)
        self.assertTrue('directed' in a_names)

        for x in ne.get_attributes():
            if x.get_name() == 'citation':
                self.assertEqual(3, len(x.get_value()))
                self.assertTrue('1' in x.get_value())
                self.assertTrue('2' in x.get_value())
                self.assertTrue('3' in x.get_value())






