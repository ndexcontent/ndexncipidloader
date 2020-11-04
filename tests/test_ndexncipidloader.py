#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `NDExNciPidLoader` class."""

import os
import tempfile
import shutil

import unittest
import mock
from mock import MagicMock
from ndex2.nice_cx_network import NiceCXNetwork

import ndexncipidloader
from ndexncipidloader.ndexloadncipid import NDExNciPidLoader
from ndexncipidloader.ndexloadncipid import NetworkIssueReport
from ndexutil.config import NDExUtilConfig
from ndexncipidloader import ndexloadncipid
from ndexncipidloader.exceptions import NDExNciPidLoaderError


class Param(object):
    """
    Dummy object
    """
    pass


class TestNDExNciPidLoader(unittest.TestCase):
    """Tests for `NDExNciPidLoader` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_parse_config(self):

        temp_dir = tempfile.mkdtemp()
        try:
            p = Param()
            p.profile = 'foo'
            conf = os.path.join(temp_dir, 'blah.conf')
            p.conf = conf
            with open(conf, 'w') as f:
                f.write('[' + p.profile + ']' + '\n')
                f.write(NDExUtilConfig.USER + ' = bob\n')
                f.write(NDExUtilConfig.PASSWORD + ' = pass\n')
                f.write(NDExUtilConfig.SERVER + ' = server\n')
                f.flush()

            loader = NDExNciPidLoader(p)
            loader._parse_config()
            self.assertEqual('bob', loader._user)
            self.assertEqual('pass', loader._pass)
            self.assertEqual('server', loader._server)
        finally:
            shutil.rmtree(temp_dir)

    def test_parse_load_plan(self):
        p = Param()
        p.loadplan = ndexloadncipid.get_load_plan()
        loader = NDExNciPidLoader(p)
        loader._parse_load_plan()
        self.assertTrue(isinstance(loader._loadplan, dict))

    def test_get_style_template(self):
        p = Param()
        p.style = ndexloadncipid.get_style()
        loader = NDExNciPidLoader(p)
        loader._load_style_template()
        self.assertTrue(loader._template is not None)

    def test_normalize_context_prefixes(self):
        loader = NDExNciPidLoader(None)
        self.assertEqual(None,
                         loader._normalize_context_prefixes(None))

        # try with value that shouldn't change
        self.assertEqual('foo',
                         loader._normalize_context_prefixes('foo'))

        # try uniprot knowledgebase:
        self.assertEqual('uniprot:hi',
                         loader._normalize_context_prefixes('uniprot '
                                                            'knowledgebase:'
                                                            'hi'))
        # try kegg compound:
        self.assertEqual('kegg.compound:hi',
                         loader._normalize_context_prefixes('kegg compound:'
                                                            'hi'))
        # try UniProt:
        self.assertEqual('uniprot:yo',
                         loader._normalize_context_prefixes('UniProt:'
                                                            'yo'))

    def test_set_wasderivedfrom(self):
        net = NiceCXNetwork()
        net.set_name('foo')
        loader = NDExNciPidLoader(None)
        loader._set_wasderivedfrom(net)
        derived_attr = ndexloadncipid.DERIVED_FROM_ATTRIB
        self.assertEqual('<a href="ftp://' +
                         ndexloadncipid.DEFAULT_FTP_HOST +
                         '/' + ndexloadncipid.DEFAULT_FTP_DIR + '/' +
                         'foo.owl.gz">foo.owl.gz</a>',
                         net.get_network_attribute(derived_attr)['v'])

    def test_set_normalizationversion(self):
        net = NiceCXNetwork()
        net.set_name('foo')
        loader = NDExNciPidLoader(None)
        loader._set_normalization_version(net)
        norm_attr = ndexloadncipid.NORMALIZATIONVERSION_ATTRIB
        self.assertEqual('0.1',
                         net.get_network_attribute(norm_attr)['v'])

    def test_set_generatedby_in_network_attributes(self):
        net = NiceCXNetwork()
        net.set_name('foo')
        loader = NDExNciPidLoader(None)
        loader._set_generatedby_in_network_attributes(net)
        norm_attr = ndexloadncipid.GENERATED_BY_ATTRIB
        self.assertTrue(' ' + str(ndexncipidloader.__version__) in
                        net.get_network_attribute(norm_attr)['v'])

    def test_set_type(self):
        loader = NDExNciPidLoader(None)
        net = NiceCXNetwork()
        net.set_name('foo')
        loader._set_type(net)
        self.assertEqual(['pathway'],
                         net.get_network_attribute('networkType')['v'])

        net.set_name(ndexloadncipid.COMPLETE_INTERACTION_NAME)
        loader._set_type(net)
        self.assertEqual(['pathway', 'interactome'],
                         net.get_network_attribute('networkType')['v'])

    def test_set_iconurl(self):
        p = Param()
        p.iconurl='hi'
        loader = NDExNciPidLoader(p)
        net = NiceCXNetwork()
        loader._set_iconurl(net)
        self.assertEqual('hi',
                         net.get_network_attribute('__iconurl')['v'])

    def test_add_node_types_in_network_to_report_empty_network(self):
        loader = NDExNciPidLoader(None)
        report = NetworkIssueReport('foo')
        net = NiceCXNetwork()
        net.set_name('foo')
        loader._add_node_types_in_network_to_report(net, report)
        self.assertEqual(set(), report.get_nodetypes())

    def test_set_network_attributes_from_style_network_description_none(self):
        net = NiceCXNetwork()
        net.set_name('foo')

        templatenet = NiceCXNetwork()
        templatenet.set_name('well')
        templatenet.set_network_attribute('organism', values='hi', type='string')

        loader = NDExNciPidLoader(None)
        loader._template = templatenet

        res = loader._set_network_attributes_from_style_network(net)
        self.assertEqual(1, len(res))
        self.assertTrue('description network' in res[0])

        self.assertEqual('', net.get_network_attribute('description')['v'])
        self.assertEqual('hi', net.get_network_attribute('organism')['v'])

    def test_set_network_attributes_from_style_network_description_set_org_not(self):
        net = NiceCXNetwork()
        net.set_name('foo')

        templatenet = NiceCXNetwork()
        templatenet.set_name('well')
        templatenet.set_network_attribute('description', values='hi', type='string')
        loader = NDExNciPidLoader(None)
        loader._template = templatenet

        res = loader._set_network_attributes_from_style_network(net)
        self.assertEqual(1, len(res))
        self.assertTrue('organism network' in res[0])

        self.assertEqual('hi', net.get_network_attribute('description')['v'])
        self.assertEqual(None, net.get_network_attribute('organism'))

    def test_set_network_attributes_from_style_network_complete_net(self):
        net = NiceCXNetwork()
        net.set_name(ndexloadncipid.COMPLETE_INTERACTION_NAME)

        templatenet = NiceCXNetwork()
        templatenet.set_name('well')
        templatenet.set_network_attribute('description', values='hi', type='string')
        templatenet.set_network_attribute('organism', values='some', type='string')
        loader = NDExNciPidLoader(None)
        loader._template = templatenet

        res = loader._set_network_attributes_from_style_network(net)
        self.assertEqual(0, len(res))

        self.assertTrue('This network' in net.get_network_attribute('description')['v'])
        self.assertEqual('some', net.get_network_attribute('organism')['v'])

    def test_update_network_system_properties_success(self):
        net = NiceCXNetwork()
        net.set_name('foo')
        loader = NDExNciPidLoader(None)
        loader._ndex = MagicMock()
        loader._ndex.set_network_system_properties = MagicMock(return_value='')
        perm_dict = {'index_level': 'ALL',
                     'showcase': True}
        loader._update_network_system_properties('http://public.ndexbio.org'
                                                 '/v2/network/someid')
        loader._ndex.set_network_system_properties.assert_called_with('someid',
                                                                      perm_dict)

    def test_update_network_system_properties_failure(self):
        net = NiceCXNetwork()
        net.set_name('foo')
        loader = NDExNciPidLoader(None)
        loader._ndex = MagicMock()
        loader._networksystemproperty_retry = 1
        loader._networksystemproperty_wait = 0
        loader._ndex.set_network_system_properties = MagicMock(return_value='error')
        perm_dict = {'index_level': 'ALL',
                     'showcase': True}
        res = loader._update_network_system_properties('http://public.'
                                                       'ndexbio.org'
                                                       '/v2/network/someid')
        self.assertTrue('After 1 retries' in res)
        loader._ndex.set_network_system_properties.assert_called_with('someid',
                                                                      perm_dict)

    def test_update_network_system_properties_exception_then_okay(self):
        net = NiceCXNetwork()
        net.set_name('foo')
        loader = NDExNciPidLoader(None)
        loader._ndex = MagicMock()
        loader._networksystemproperty_wait = 0
        loader._ndex.set_network_system_properties = MagicMock(side_effect=[Exception('foo'), ''])
        perm_dict = {'index_level': 'ALL',
                     'showcase': True}
        res = loader._update_network_system_properties('http://public.'
                                                       'ndexbio.org'
                                                       '/v2/network/someid')
        self.assertEqual(None, res)
        loader._ndex.set_network_system_properties.assert_called_with('someid',
                                                                      perm_dict)

        cnt = loader._ndex.set_network_system_properties.call_count
        self.assertEqual(2, cnt)

    def test_apply_simple_spring_layout(self):
        net = NiceCXNetwork()
        for x in range(2):
            net.create_node('node' + str(x))
        loader = NDExNciPidLoader(None)
        loader._apply_simple_spring_layout(net)
        res = net.get_opaque_aspect('cartesianLayout')
        self.assertEqual(2, len(res))
        self.assertTrue('node' in res[0])
        self.assertTrue(res[0]['node'] >= 0)
        self.assertTrue('x' in res[0])
        self.assertTrue('y' in res[0])
        self.assertTrue('node' in res[1])
        self.assertTrue(res[1]['node'] >= 1)
        self.assertTrue('x' in res[1])
        self.assertTrue('y' in res[1])

        net = NiceCXNetwork()
        for x in range(10):
            net.create_node('node' + str(x))
        loader._apply_simple_spring_layout(net)
        self.assertEqual(10,
                         len(net.get_opaque_aspect('cartesianLayout')))

        net = NiceCXNetwork()
        for x in range(20):
            net.create_node('node' + str(x))
        loader._apply_simple_spring_layout(net)
        self.assertEqual(20,
                         len(net.get_opaque_aspect('cartesianLayout')))

        net = NiceCXNetwork()
        for x in range(100):
            net.create_node('node' + str(x))
        loader._apply_simple_spring_layout(net)
        self.assertEqual(100,
                         len(net.get_opaque_aspect('cartesianLayout')))

    def test_apply_cytoscape_layout_ping_failed(self):
        p = MagicMock()
        p.layout = 'grid'
        mockpy4 = MagicMock()
        mockpy4.cytoscape_ping = MagicMock(side_effect=Exception('error'))
        loader = NDExNciPidLoader(p, py4cyto=mockpy4)
        net = NiceCXNetwork()
        try:
            loader._apply_cytoscape_layout(net)
            self.fail('Expected NDExNciPidLoaderError')
        except NDExNciPidLoaderError as e:
            self.assertEqual('Cytoscape needs to be running '
                             'to run layout: grid', str(e))

    def test_apply_cytoscape_layout_networks_not_in_dict(self):
        temp_dir = tempfile.mkdtemp()
        try:
            p = MagicMock()
            p.layout = 'grid'
            p.sifdir = temp_dir
            mockpy4 = MagicMock()
            mockpy4.import_network_from_file = MagicMock(return_value={})
            loader = NDExNciPidLoader(p, py4cyto=mockpy4)
            net = NiceCXNetwork()
            for x in range(10):
                net.create_node('node' + str(x))
            try:
                loader._apply_cytoscape_layout(net)
                self.fail('Expected NDExNciPidLoaderError')
            except NDExNciPidLoaderError as e:
                self.assertTrue(str(e).startswith('Error network view'))
        finally:
            shutil.rmtree(temp_dir)

    def test_apply_cytoscape_layout_networks_success(self):
        temp_dir = tempfile.mkdtemp()
        try:
            p = MagicMock()
            p.layout = 'grid'
            p.sifdir = temp_dir
            mockpy4 = MagicMock()
            imp_res = {'networks': ['netid']}
            mockpy4.import_network_from_file = MagicMock(return_value=imp_res)
            mockpy4.export_network = MagicMock(return_value='')
            loader = NDExNciPidLoader(p, py4cyto=mockpy4)
            loader._ndexextra.extract_layout_aspect_from_cx = MagicMock(return_value={'cartesianLayout': []})
            net = NiceCXNetwork()
            for x in range(10):
                net.create_node('node' + str(x))
            loader._apply_cytoscape_layout(net)
            self.assertEqual([{'cartesianLayout': []}],
                             net.get_opaque_aspect('cartesianLayout'))
        finally:
            shutil.rmtree(temp_dir)



