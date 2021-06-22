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
from ndexindraloader.indra import Indra


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

    def test_get_lists_of_edges_to_merge(self):
        annotator = INDRAAnnotator()
        net = NiceCXNetwork()

        node_one = net.create_node('node1')
        node_two = net.create_node('node2')
        node_three = net.create_node('node3')

        edge_one = net.create_edge(edge_source=node_one,
                                   edge_target=node_two)
        net.set_edge_attribute(edge=edge_one, attribute_name=Indra.SOURCE,
                               values='NCI PID')
        edge_two = net.create_edge(edge_source=node_one,
                                   edge_target=node_two)

        edge_three = net.create_edge(edge_source=node_two,
                                     edge_target=node_three)

        edge_four = net.create_edge(edge_source=node_one,
                                     edge_target=node_three)
        net.set_edge_attribute(edge=edge_four, attribute_name=Indra.SOURCE,
                               values='NCI PID')

        edge_five = net.create_edge(edge_source=node_three,
                                     edge_target=node_one)

        nci_edges = annotator._get_nci_pid_edge_objs(network=net)
        self.assertEqual(2, len(nci_edges))
        res = annotator._get_lists_of_edges_to_merge(net_cx=net,
                                                     nci_pid_edges=nci_edges)
        self.assertEqual(2, len(res))
        self.assertEqual(2, len(res[0]))
        self.assertEqual(2, len(res[1]))

        # TODO finish testing this method

    def test_split_edges_into_indra_and_nci_pid(self):
        annotator = INDRAAnnotator()

        edges = [NetworkEdge(edge_id=0, interaction=None),
                 NetworkEdge(edge_id=1, interaction='foo'),
                 NetworkEdge(edge_id=2, interaction='interacts with')]

        res = annotator._split_edges_into_indra_and_nci_pid(edges=edges)
        self.assertEqual(2, res[0].get_id())
        self.assertEqual(2, len(res[1]))
        self.assertTrue(res[1][0].get_id() != 2)
        self.assertTrue(res[1][1].get_id() != 2)



    def test_is_edge_directed(self):
        annotator = INDRAAnnotator()

        # try where edge has no attributes
        edge = NetworkEdge(edge_id=0, source_node_id=1,
                           target_node_id=2)

        self.assertFalse(annotator._is_edge_directed(edge=edge))

        # try where edge does not have directed attribute
        edge.set_attributes([Attribute(name='foo', value='blah'),
                             Attribute(name='foo2', value='blah2')])
        self.assertFalse(annotator._is_edge_directed(edge=edge))

        # try where directed attribute is False
        edge.set_attributes([Attribute(name='foo', value='blah'),
                             Attribute(name=Indra.DIRECTED, value=False,
                                       data_type='boolean')])
        self.assertFalse(annotator._is_edge_directed(edge=edge))

        # try where directed attribute is True
        edge.set_attributes([Attribute(name='foo', value='blah'),
                             Attribute(name=Indra.DIRECTED, value=True,
                                       data_type='boolean')])
        self.assertTrue(annotator._is_edge_directed(edge=edge))

    def test_get_nci_citations(self):
        annotator = INDRAAnnotator()

        edge = NetworkEdge(edge_id=0, source_node_id=1,
                           target_node_id=2)

        # try where no citation attribute exists
        self.assertEqual([], annotator._get_nci_citations(nci_edge=edge))

        edge.set_attributes([Attribute(name='foo', value='blah'),
                             Attribute(name='foo2', value='blah2')])
        self.assertEqual([], annotator._get_nci_citations(nci_edge=edge))

        # try where citation does exist
        edge.get_attributes().append(Attribute(name='citation',
                                               value=['1', '2'],
                                               data_type='list_of_string'))
        self.assertEqual(['1', '2'],
                         annotator._get_nci_citations(nci_edge=edge))

    def test_update_interaction_attribute_on_empty_edge(self):
        annotator = INDRAAnnotator()

        # try where no interactions exist and nor do directed or
        # reverse directed and both directed, reverse are true
        edge = NetworkEdge(edge_id=1, source_node_id=2, target_node_id=3)
        res = annotator._update_interaction_attribute(edge=edge,
                                                      name='SOURCE => TARGET',
                                                      interaction='breaks',
                                                      citations=['pubmed:1',
                                                                 'pubmed:2'],
                                                      update_directed=True,
                                                      update_reversedirected=True)

        self.assertEqual([], res)

        # check directed
        dir_attr = edge.get_attribute_by_name(Indra.DIRECTED)
        self.assertEqual(True, dir_attr.get_value())
        self.assertEqual('boolean', dir_attr.get_data_type())

        # check reverse directed
        rdir_attr = edge.get_attribute_by_name(Indra.REVERSE_DIRECTED)
        self.assertEqual(True, rdir_attr.get_value())
        self.assertEqual('boolean', rdir_attr.get_data_type())

        # check interaction
        int_attr = edge.get_attribute_by_name('SOURCE => TARGET')
        self.assertEqual('list_of_string', int_attr.get_data_type())
        self.assertEqual(['breaks(NCI PID - pubmed:1,pubmed:2)'],
                         int_attr.get_value())

    def test_update_interaction_attribute_on_undirected_edge(self):
        annotator = INDRAAnnotator()

        # try where no interactions exist and nor do directed or
        # reverse directed and both directed, reverse are true
        edge = NetworkEdge(edge_id=1, source_node_id=2, target_node_id=3)
        edge.set_attributes([Attribute(name=Indra.DIRECTED, value=False,
                                      data_type='boolean'),
                             Attribute(name='SOURCE - TARGET',
                                       value=['blah(1)'],
                                       data_type='list_of_string'),
                             Attribute(name=Indra.REVERSE_DIRECTED,
                                       value=False, data_type='boolean')])
        res = annotator._update_interaction_attribute(edge=edge,
                                                      name='SOURCE => TARGET',
                                                      interaction='breaks',
                                                      citations=['pubmed:1',
                                                                 'pubmed:2'],
                                                      update_directed=True,
                                                      update_reversedirected=True)

        self.assertEqual([], res)

        # check directed
        dir_attr = edge.get_attribute_by_name(Indra.DIRECTED)
        self.assertEqual(True, dir_attr.get_value())
        self.assertEqual('boolean', dir_attr.get_data_type())

        # check reverse directed
        rdir_attr = edge.get_attribute_by_name(Indra.REVERSE_DIRECTED)
        self.assertEqual(True, rdir_attr.get_value())
        self.assertEqual('boolean', rdir_attr.get_data_type())

        # check interaction
        int_attr = edge.get_attribute_by_name('SOURCE => TARGET')
        self.assertEqual('list_of_string', int_attr.get_data_type())
        self.assertEqual(['breaks(NCI PID - pubmed:1,pubmed:2)'],
                         int_attr.get_value())

        int_attr = edge.get_attribute_by_name('SOURCE - TARGET')
        self.assertEqual('list_of_string', int_attr.get_data_type())
        self.assertEqual(['blah(1)'],
                         int_attr.get_value())

    def test_get_nci_pid_edge_ids(self):
        annotator = INDRAAnnotator()
        net = NiceCXNetwork()

        node_one = net.create_node('node1')
        node_two = net.create_node('node2')
        node_three = net.create_node('node3')

        edge_one = net.create_edge(edge_source=node_one,
                                   edge_target=node_two)
        net.set_edge_attribute(edge=edge_one,
                               attribute_name=Indra.SOURCE,
                               values='NCI PID')
        edge_two = net.create_edge(edge_source=node_one,
                                   edge_target=node_three)
        net.set_edge_attribute(edge=edge_two,
                               attribute_name=Indra.SOURCE,
                               values='NCI PID')

        edge_three = net.create_edge(edge_source=node_two,
                                     edge_target=node_three)

        edge_ids = annotator._get_nci_pid_edge_objs(network=net)

        self.assertEqual(2, len(edge_ids))

        e_dict = {}
        for e_tup in edge_ids:
            e_dict[e_tup[0]] = e_tup[1]
        self.assertTrue(edge_one in e_dict)
        self.assertEqual({'@id': 0, 's': 0, 't': 1}, e_dict[edge_one])
        self.assertTrue(edge_two in e_dict)
        self.assertEqual({'@id': 1, 's': 0, 't': 2}, e_dict[edge_two])

