#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ProteinFamilyNodeMemberRemover` class."""

import os
import logging
import unittest
from ndex2.nice_cx_network import NiceCXNetwork
import ndex2

from ndexncipidloader.ndexloadncipid import ProteinFamilyNodeMemberRemover
from ndexncipidloader.network import NetworkEdgeFactory
from ndexncipidloader.network import NetworkEdge
from ndexncipidloader.network import Attribute


class TestProteinFamilyNodeMemberRemover(unittest.TestCase):
    """Tests for `ProteinFamilyNodeMemberRemover` class."""

    def get_bmp_receptor_signaling_file(self):
        """
        Gets path to bmp receptor signaling file
        :return:
        """
        return os.path.join(os.path.dirname(__file__),
                            'data', 'bmp_receptor_signaling.cx')

    def setUp(self):
        """Set up test fixtures, if any."""
        # set the logging level to debug cause some code
        # will only run in this logging level
        x = logging.getLogger('ndexncipidloader.ndexloadncipid')
        x.setLevel(logging.DEBUG)

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_description(self):
        remover = ProteinFamilyNodeMemberRemover()
        self.assertTrue('Removes nodes that '
                        'are part' in remover.get_description())

    def test_get_edge_key(self):
        ne = NetworkEdge(source_node_id=1, target_node_id=2,
                         source_node_name='src', target_node_name='target',
                         interaction='activates')
        remover = ProteinFamilyNodeMemberRemover()

        # try with no attributes
        e_key, d_val = remover._get_edge_key(ne)
        self.assertEqual('s=1, t=2, i=activates, directed=False', e_key)
        self.assertEqual(False, d_val)
        e_key, d_val = remover._get_edge_key(ne, flip_src_target=True)
        self.assertEqual('s=2, t=1, i=activates, directed=False', e_key)
        self.assertEqual(False, d_val)

        # try with attributes but no match
        ne.set_attributes([Attribute(name='foo', value='1')])
        e_key, d_val = remover._get_edge_key(ne)
        self.assertEqual('s=1, t=2, i=activates, directed=False', e_key)
        self.assertEqual(False, d_val)
        e_key, d_val = remover._get_edge_key(ne, flip_src_target=True)
        self.assertEqual('s=2, t=1, i=activates, directed=False', e_key)
        self.assertEqual(False, d_val)

        # try with directed attribute True
        ne.set_attributes([Attribute(name='directed', value=True,
                                     data_type='boolean')])
        e_key, d_val = remover._get_edge_key(ne)
        self.assertEqual('s=1, t=2, i=activates, directed=True', e_key)
        self.assertEqual(True, d_val)
        e_key, d_val = remover._get_edge_key(ne, flip_src_target=True)
        self.assertEqual('s=2, t=1, i=activates, directed=True', e_key)
        self.assertEqual(True, d_val)

        # try with directed attribute False
        ne.set_attributes([Attribute(name='directed', value=False,
                                     data_type='boolean')])
        e_key, d_val = remover._get_edge_key(ne)
        self.assertEqual('s=1, t=2, i=activates, directed=False', e_key)
        self.assertEqual(False, d_val)
        e_key, d_val = remover._get_edge_key(ne, flip_src_target=True)
        self.assertEqual('s=2, t=1, i=activates, directed=False', e_key)
        self.assertEqual(False, d_val)

    def test_get_edge_key(self):
        remover = ProteinFamilyNodeMemberRemover()

        # try with 1 edge no attributes
        ne = NetworkEdge(source_node_id=3, target_node_id=2,
                         interaction='activates')
        res = remover._get_edge_dict([ne])
        self.assertEqual(1, len(res))
        key = 's=3, t=2, i=activates, directed=False'
        self.assertTrue(key in res)
        self.assertEqual(1, len(res[key]))
        self.assertTrue(ne.get_source_node_id() ==
                        res[key][0].get_source_node_id())

        # try 3 edges where two edges have source/target are flipped and
        # directed is same value
        ne_one = NetworkEdge(source_node_id=1, target_node_id=2,
                             interaction='activates',
                             attributes=[Attribute(name='directed',
                                                   value=False,
                                                   data_type='boolean')])
        ne_two = NetworkEdge(source_node_id=2, target_node_id=1,
                             interaction='activates',
                             attributes=[Attribute(name='directed',
                                                   value=False,
                                                   data_type='boolean')])
        res = remover._get_edge_dict([ne_one, ne_two, ne])

        self.assertEqual(2, len(res))
        self.assertTrue(key in res)
        self.assertEqual(1, len(res[key]))
        self.assertTrue(ne.get_source_node_id() ==
                        res[key][0].get_source_node_id())
        mkey = 's=1, t=2, i=activates, directed=False'
        self.assertTrue(mkey in res)
        self.assertEqual(2, len(res[mkey]))





    def test_update_no_family_members(self):
        net = NiceCXNetwork()

        node1 = net.create_node('node1')
        node2 = net.create_node('node2')
        net.create_edge(edge_source=node1, edge_target=node2,
                        edge_interaction='activates')
        remover = ProteinFamilyNodeMemberRemover()
        remover.update(net)

        self.assertEqual(2, len(net.get_nodes()))
        self.assertEqual(1, len(net.get_edges()))

    def test_update_with_bmp_receptor_signaling(self):
        net = ndex2.\
            create_nice_cx_from_file(self.get_bmp_receptor_signaling_file())

        self.assertEqual(43, len(net.get_nodes()))
        self.assertEqual(224, len(net.get_edges()))
        remover = ProteinFamilyNodeMemberRemover()
        remover.update(net)

        self.assertEqual(39, len(net.get_nodes()))
        self.assertEqual(168, len(net.get_edges()))

        # verify TAB1, TAB2, RGMB, & RGMA nodes have been removed
        node_name_to_id = remover._get_node_name_to_id_dict(network=net)
        self.assertTrue('RGMB' not in node_name_to_id)
        self.assertTrue('RGMA' not in node_name_to_id)
        self.assertTrue('TAB1' not in node_name_to_id)
        self.assertTrue('TAB2' not in node_name_to_id)

        nef = NetworkEdgeFactory()
        rgm_family = node_name_to_id['RGM family']
        edges = nef.get_all_edges_connected_to_node(net_cx=net,
                                                    node_id=rgm_family)
        self.assertEqual(16, len(edges))

        tab_family = node_name_to_id['TAB family']
        edges = nef.get_all_edges_connected_to_node(net_cx=net,
                                                    node_id=tab_family)
        self.assertEqual(10, len(edges))






