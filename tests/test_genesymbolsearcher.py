#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `GeneSymbolSearcher` class."""

import os
import tempfile
import shutil

import unittest
import mock
from mock import MagicMock

from ndexncipidloader.loadndexncipidloader import GeneSymbolSearcher


class TestGeneSymbolSearcher(unittest.TestCase):
    """Tests for `GeneSymbolSearcher` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_symbol_none_passed_in(self):
        searcher = GeneSymbolSearcher()
        self.assertEqual(None, searcher.get_symbol(None))

    def test_symbol_not_in_cache_no_hit(self):
        mock = MagicMock()
        mock.query = MagicMock(return_value={'max_score': None,
                                             'took': 5, 'total': 0,
                                             'hits': []})
        searcher = GeneSymbolSearcher(bclient=mock)
        self.assertEqual(None, searcher.get_symbol('haha'))

        mock.query.assert_called_with('haha')

    def test_symbol_not_in_cache_no_hit_and_none_resturned(self):
        mock = MagicMock()
        mock.query = MagicMock(return_value=None)
        searcher = GeneSymbolSearcher(bclient=mock)
        self.assertEqual(None, searcher.get_symbol('haha'))

        mock.query.assert_called_with('haha')

    def test_symbol_not_in_cache_no_hit_and_hit_name_is_none(self):
        mock = MagicMock()
        mock.query = MagicMock(return_value={'total': 1,
                                             'hits': [{'symbol': None}]})
        searcher = GeneSymbolSearcher(bclient=mock)
        self.assertEqual(None, searcher.get_symbol('haha'))

        mock.query.assert_called_with('haha')

    def test_symbol_not_in_cache_no_hit_total_incorrect(self):
        mock = MagicMock()
        mock.query = MagicMock(return_value={'total': 1,
                                             'hits': []})
        searcher = GeneSymbolSearcher(bclient=mock)
        self.assertEqual(None, searcher.get_symbol('haha'))

        mock.query.assert_called_with('haha')

    def test_symbol_not_in_cache_no_but_got_hit_with_cache_check(self):
        mock = MagicMock()
        mock.query = MagicMock(return_value={'max_score': 437.58682, 'took': 62,
                                             'total': 5031,
                                             'hits': [{'_id': '7157',
                                                       '_score': 437.58682,
                                                       'entrezgene': '7157',
                                                       'name': 'tumor protein p53',
                                                       'symbol': 'TP53',
                                                       'taxid': 9606},
                                                      {'_id': '24842',
                                                       '_score': 306.22318,
                                                       'entrezgene': '24842',
                                                       'name': 'tumor protein p53',
                                                       'symbol': 'Tp53',
                                                       'taxid': 10116},
                                                      {'_id': '109394672',
                                                       '_score': 296.02454,
                                                       'entrezgene': '109394672',
                                                       'name': 'tumor protein p53',
                                                       'symbol': 'TP53',
                                                       'taxid': 186990},
                                                      {'_id': '113633022',
                                                       '_score': 296.02454,
                                                       'entrezgene': '113633022',
                                                       'name': 'tumor protein p53',
                                                       'symbol': 'TP53',
                                                       'taxid': 90247},
                                                      {'_id': '102169621', '_score': 296.02454, 'entrezgene': '102169621', 'name': 'tumor protein p53', 'symbol': 'TP53', 'taxid': 9925}, {'_id': '113878373', '_score': 296.02454, 'entrezgene': '113878373', 'name': 'tumor protein p53', 'symbol': 'TP53', 'taxid': 30522}, {'_id': '101285670', '_score': 296.02454, 'entrezgene': '101285670', 'name': 'tumor protein p53', 'symbol': 'TP53', 'taxid': 9733}, {'_id': '105819395', '_score': 296.02454, 'entrezgene': '105819395', 'name': 'tumor protein p53', 'symbol': 'TP53', 'taxid': 379532}, {'_id': 'ENSMMMG00000019747', '_score': 296.02454, 'name': 'tumor protein p53', 'symbol': 'TP53', 'taxid': 9994}, {'_id': '100583326', '_score': 296.02454, 'entrezgene': '100583326', 'name': 'tumor protein p53', 'symbol': 'TP53', 'taxid': 61853}]})
        searcher = GeneSymbolSearcher(bclient=mock)
        self.assertEqual('TP53', searcher.get_symbol('tp53'))
        self.assertEqual('TP53', searcher.get_symbol('tp53'))
        mock.query.assert_called_once_with('tp53')
