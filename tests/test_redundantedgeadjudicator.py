#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `RedundantEdgeAdjudicator` class."""

import os
import tempfile
import shutil

import unittest
import mock
from mock import MagicMock

from ndexncipidloader.loadndexncipidloader import RedundantEdgeAdjudicator
from ndex2.nice_cx_network import NiceCXNetwork

class TestRedundantEdgeAdjudicator(unittest.TestCase):
    """Tests for `RedundantEdgeAdjudicator` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_remove_nonexistant_edge(self):
        adjud = RedundantEdgeAdjudicator()
        net = NiceCXNetwork()
        adjud._remove_edge(net, 1)

    def test_remove_edge_no_attributes(self):
        adjud = RedundantEdgeAdjudicator()
        net = NiceCXNetwork()
        edgeid = net.create_edge(edge_source=0, edge_target=1,
                                 edge_interaction='foo')
        self.assertEqual('foo', net.get_edge(edgeid)['i'])
        adjud._remove_edge(net, edgeid)
        self.assertEqual(None, net.get_edge(edgeid))

    def test_remove_edge_with_attributes(self):
        adjud = RedundantEdgeAdjudicator()
        net = NiceCXNetwork()
        edgeid = net.create_edge(edge_source=0, edge_target=1,
                                 edge_interaction='foo')
        net.set_edge_attribute(edgeid, 'attr1', 'someval')
        self.assertEqual('someval', net.get_edge_attribute(edgeid,
                                                           'attr1')['v'])
        self.assertEqual('foo', net.get_edge(edgeid)['i'])
        adjud._remove_edge(net, edgeid)
        self.assertEqual(None, net.get_edge(edgeid))
        self.assertEqual((None, None),
                         net.get_edge_attribute(edgeid, 'attr1'))
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
