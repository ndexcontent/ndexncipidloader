#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `DirectedEdgeSetter` class."""

import os
import tempfile
import shutil

import unittest
import mock
from mock import MagicMock

from ndexncipidloader import loadndexncipidloader
from ndexncipidloader.loadndexncipidloader import DirectedEdgeSetter
from ndex2.nice_cx_network import NiceCXNetwork

class TestDirectedEdgeSetter(unittest.TestCase):
    """Tests for `DirectedEdgeSetter` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_with_none_as_network(self):
        setter = DirectedEdgeSetter()
        self.assertEqual(['Network is None'],
                         setter.update_edge_direction(None))

    def test_with_empty_network(self):
        setter = DirectedEdgeSetter()
        net = NiceCXNetwork()
        self.assertEqual([], setter.update_edge_direction(net))

    def test_with_network_with_edge_not_supposed_to_be_directed(self):
        setter = DirectedEdgeSetter()
        net = NiceCXNetwork()
        edgeid = net.create_edge(edge_source=0,
                                 edge_target=1,
                                 edge_interaction='foo')
        self.assertEqual([], setter.update_edge_direction(net))
        directed = DirectedEdgeSetter.DIRECTED_ATTRIB
        self.assertEqual(False,
                         net.get_edge_attribute(edgeid,
                                                directed)['v'])

    def test_with_network_with_edges_that_are_directed(self):
        setter = DirectedEdgeSetter()
        net = NiceCXNetwork()
        directed = DirectedEdgeSetter.DIRECTED_ATTRIB
        cntr = 1
        for entry in loadndexncipidloader.DIRECTED_INTERACTIONS:

            edgeid = net.create_edge(edge_source=0,
                                     edge_target=cntr,
                                     edge_interaction=entry)
            self.assertEqual([], setter.update_edge_direction(net))
            self.assertEqual(True,
                             net.get_edge_attribute(edgeid,
                                                    directed)['v'])
            cntr = cntr + 1

"""
    def test_symbol_not_in_cache_no_hit(self):
        mock = MagicMock()
        mock.query = MagicMock(return_value={'max_score': None,
                                             'took': 5, 'total': 0,
                                             'hits': []})
        searcher = GeneSymbolSearcher(bclient=mock)
        self.assertEqual(None, searcher.get_symbol('haha'))

        mock.query.assert_called_with('haha')
"""
