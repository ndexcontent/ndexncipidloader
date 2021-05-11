# -*- coding: utf-8 -*-

import logging


logger = logging.getLogger(__name__)


class Attribute(object):
    """
    Object that represents an Attribute
    """
    def __init__(self, name=None, value=None, data_type=None):
        """
        Constructor
        """
        self._name = name
        self._value = value
        self._data_type = data_type

    def get_name(self):
        """
        Gets edge annotation/attribute name

        :return:
        """
        return self._name

    def set_name(self, name):
        """
        Sets edge annotation/attribute name

        :param name:
        :return:
        """
        self._name = name

    def get_value(self):
        """
        Gets edge annotation/attribute value

        :return:
        """
        return self._value

    def set_value(self, value):
        """
        Sets edge annotation/attribute name

        :param value:
        :return:
        """
        self._value = value

    def get_data_type(self):
        """
        Gets edge annotation/attribute data type

        :return:
        """
        return self._data_type

    def set_data_type(self, data_type):
        """
        Sets edge annotation/attribute data type

        :param data_type:
        :return:
        """
        self._data_type = data_type


class NetworkEdgeAttribute(Attribute):
    """
    Object that represents an edge annotation
    """
    def __init__(self, name=None, value=None, data_type=None):
        """
        Constructor
        """
        super(NetworkEdgeAttribute, self).__init__(name=name,
                                                   value=value,
                                                   data_type=data_type)
        self._equalityfailreason = None

    def add_attribute_to_edge(self, net_cx=None, edge_id=None):
        """
        Adds this edge attribute

        :param net_cx: network to add edge attribute to
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param edge_id: id of edge on **net_cx** network
        :type edge_id: int
        :return:
        """
        net_cx.set_edge_attribute(edge_id, self._name,
                                  self._value,
                                  type=self._data_type)

    def __str__(self):
        """
        Returns str representation of object
        :return:
        """
        return 'name=' + str(self._name) + ', values=' +\
               str(self._value) + ', type=' +\
               str(self._data_type)

    def __eq__(self, other):
        """
        Compares self with **other** for equality
        :param other:
        :type other: :py:class:`NetworkEdgeAttribute`
        :return: True if match False otherwise
        """
        if other is None:
            self._equalityfailreason = 'Other is None'
            return False
        if self._name != other.get_name():
            self._equalityfailreason = str(self._name) +\
                                       ' name does not match ' +\
                                       str(other.get_name())
            return False

        if self._value != other.get_value():
            if isinstance(self._value, list) and\
                 isinstance(other.get_value(), list):
                my_values = sorted(list(set(self._value)))
                other_values = sorted(list(set(other.get_value())))
                if my_values != other_values:
                    self._equalityfailreason = str(self._value) + \
                                               ' value does not match ' + \
                                               str(other.get_value())
                    return False
            else:
                self._equalityfailreason = str(self._value) +\
                                           ' value does not match ' +\
                                           str(other.get_value())
                return False
        if self._data_type != other.get_data_type():
            if self._data_type == 'string' and other.get_data_type() is None:
                self._equalityfailreason = None
                return True

            if self._data_type is None and other.get_data_type() == 'string':
                self._equalityfailreason = None
                return True
            self._equalityfailreason = 'data types differ'
            return False
        self._equalityfailreason = None
        return True

    def get_equality_fail_reason(self):
        """

        :return:
        """
        return self._equalityfailreason


class NetworkNodeAttribute(Attribute):
    """
    Object that represents an node annotation
    """
    def __init__(self, name=None, value=None, data_type=None):
        """
        Constructor
        """
        super(NetworkNodeAttribute, self).__init__(name=name,
                                                   value=value,
                                                   data_type=data_type)

    def add_attribute_to_node(self, net_cx=None, node_id=None):
        """
        Adds this node attribute

        :param net_cx: network to add edge attribute to
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param edge_id: id of edge on **net_cx** network
        :type edge_id: int
        :return:
        """
        net_cx.set_node_attribute(node_id, self._name,
                                  self._value,
                                  type=self._data_type)


class NetworkNode(object):
    """
    Object that represents a node
    """
    def __init__(self, node_id=None, name=None,
                 represents=None, attributes=None,
                 network_edge_factory=None):
        self._node_id = node_id
        self._name = name
        self._represents = represents
        self._attributes = attributes
        if network_edge_factory is None:
            self._nef = NetworkEdgeFactory()
        else:
            self._nef = network_edge_factory

    def get_id(self):
        """
        Gets id of node

        :return:
        :rtype: int
        """
        return self._node_id

    def set_id(self, node_id):
        """
        Sets node id

        :param node_id:
        :return:
        """
        self._node_id = node_id

    def set_name(self, name):
        """
        Sets name

        :param name:
        :return:
        """
        self._name = name

    def get_name(self):
        """
        Gets name

        :return:
        :rtype: str
        """
        return self._name

    def get_represents(self):
        """
        Gets represents
        :return:
        """
        return self._represents

    def set_represents(self, represents):
        """
        Sets represents
        :param represents:
        :return:
        """
        self._represents = represents

    def get_attributes(self):
        """
        Gets attributes
        :return:
        """
        return self._attributes

    def set_attributes(self, attributes):
        """
        Sets attributes

        :param attributes:
        :return:
        """
        self._attributes = attributes

    def remove_node_from_network(self, net_cx=None):
        """
        Removes node and any attributes on that node
        from network as well as any edges and their
        attributes connected to that node
        :param net_cx:
        :return:
        """
        edges = self._nef.get_all_edges_connected_to_node(net_cx=net_cx,
                                                          node_id=self._node_id)
        edge_cntr = 0
        if edges is not None:
            for edge in edges:
                edge.remove_edge_from_network(net_cx=net_cx)
                edge_cntr += 1
            logger.debug('Removed ' + str(edge_cntr) +
                         ' edges linked to node: ' + str(self))

        net_cx.remove_node(self._node_id)
        node_attrs = net_cx.get_node_attributes(self._node_id)
        if node_attrs is None:
            return
        attr_names = set()
        for node_attr in node_attrs:
            attr_names.add(node_attr['n'])
        for name in attr_names:
            logger.debug('Removing node attribute: ' + str(name))
            net_cx.remove_node_attribute(self._node_id, name)
        logger.debug('Removed node: ' + str(self))

    def __str__(self):
        """
        Gets string representation of node

        :return:
        """
        return '@id=' + str(self._node_id) +\
               ', name=' + str(self._name) +\
               ', represents=' + str(self._represents)


class NetworkEdge(object):
    """
    Object that represents an edge
    """
    def __init__(self, edge_id=None, source_node_id=None,
                 source_node_name=None,
                 target_node_id=None,
                 target_node_name=None, interaction=None,
                 attributes=None):
        """
        Constructor
        """
        self._edge_id = edge_id
        self._source_node_id = source_node_id
        self._target_node_id = target_node_id
        self._interaction = interaction
        self._source_node_name = source_node_name
        self._target_node_name = target_node_name
        self._attributes = attributes

    def get_id(self):
        """
        Gets id of edge

        :return:
        :rtype: int
        """
        return self._edge_id

    def set_id(self, edge_id):
        """
        Sets edge id

        :param edge_id:
        :return:
        """
        self._edge_id = edge_id

    def get_source_node_id(self):
        """
        Gets source node id

        :return:
        :rtype: int
        """
        return self._source_node_id

    def set_source_node_id(self, source_node_id):
        """
        Sets source node id

        :param source_node_id:
        :return: None
        """
        self._source_node_id = source_node_id

    def get_target_node_id(self):
        """
        Gets target node id

        :return: target node id
        :rtype: int
        """
        return self._target_node_id

    def set_target_node_id(self, target_node_id):
        """
        Sets target node id

        :param target_node_id:
        :return: None
        """
        self._target_node_id = target_node_id

    def get_source_node_name(self):
        """
        Gets source node name

        :return:
        """
        return self._source_node_name

    def set_source_node_name(self, node_name):
        """
        Sets source node name
        :param node_name:
        :type node_name: str
        :return:
        """
        self._source_node_name = node_name

    def get_target_node_name(self):
        """
        Gets target node name

        :return:
        """
        return self._target_node_name

    def set_target_node_name(self, node_name):
        """
        Sets target node name
        :param node_name:
        :type node_name: str
        :return:
        """
        self._target_node_name = node_name

    def get_interaction(self):
        """
        Gets interaction
        :return:
        """
        return self._interaction

    def set_interaction(self, interaction):
        """
        Sets interaction

        :param interaction:
        :return:
        """
        self._interaction = interaction

    def get_attributes(self):
        """
        Gets edge annotations

        :return: :py:class:`NetworkEdgeAttribute` as a list
        :rtype: list
        """
        return self._attributes

    def set_attributes(self, attributes):
        """
        Sets edge annotations
        :param attributes: list of :py:class:`NetworkEdgeAttribute`
        :type attributes: list
        :return:
        """
        self._attributes = attributes

    def add_edge_to_network(self, net_cx=None, source_node_id=None,
                            target_node_id=None):
        """
        Adds this edge and its annotations to network
        :param net_cx:
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return:
        """
        new_edge_id = net_cx.create_edge(edge_source=source_node_id,
                                         edge_target=target_node_id,
                                         edge_interaction=self._interaction)
        if self._attributes is not None:
            for annot in self._attributes:
                annot.add_attribute_to_edge(net_cx=net_cx, edge_id=new_edge_id)
        return new_edge_id

    def remove_edge_from_network(self, net_cx=None):
        """
        Removes this edge from network passed in

        :param net_cx:
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :return:
        """
        net_cx.remove_edge(self._edge_id)
        # remove edge attributes for deleted edge
        net_attrs = net_cx.get_edge_attributes(self._edge_id)
        if net_attrs is None:
            return
        attr_names = set()
        for e_obj in net_attrs:
            attr_names.add(e_obj['n'])
        attr_cntr = 0
        for name in attr_names:
            attr_cntr += 1
            net_cx.remove_edge_attribute(self._edge_id,
                                         name)
        logger.debug('Removed ' + str(attr_cntr) +
                     ' edge attributes for ' + str(self))

    def __str__(self):
        """
        Gets string representation
        :return:
        """
        return 'po=' + str(self._edge_id) + ', s=' +\
               str(self._source_node_id) + ' (' +\
               str(self._source_node_name) + '), t=' +\
               str(self._target_node_id) + ' (' +\
               str(self._target_node_name) + '), i=' +\
               str(self._interaction)


class NetworkNodeFactory(object):
    """
    Factory to create NetworkNode objects
    """
    def __init__(self):
        """
        Constructor
        """
        pass

    def get_network_node_from_network(self, net_cx=None, node_id=None):
        """
        Gets a :py:class:`NetworkNode` from
        :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param net_cx: Network to extract node from
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param node_id: id of node in `net_cx` network
        :type node_id: int
        :return: Node or `None` if not found
        :rtype: :py:class:`NetworkNode`
        """
        node_obj = net_cx.get_node(node_id)
        if node_obj is None or node_obj is (None, None):
            return None

        attributes = []
        node_attrs = net_cx.get_node_attributes(node_id)
        if node_attrs is not None:
            for n_attr in node_attrs:
                if 'd' in n_attr:
                    data_type = n_attr['d']
                else:
                    data_type = None
                n_annot = NetworkNodeAttribute(name=n_attr['n'],
                                               value=n_attr['v'],
                                               data_type=data_type)
                attributes.append(n_annot)

        if 'r' in node_obj:
            represents = node_obj['r']
        else:
            represents = None

        return NetworkNode(node_id=node_obj['@id'],
                           name=node_obj['n'],
                           represents=represents,
                           attributes=attributes)


class NetworkEdgeFactory(object):
    """
    Factory to create NetworkEdge objects
    """
    def __init__(self):
        """
        Constructor
        """
        pass

    def get_network_edge_from_network(self, net_cx=None, edge_id=None):
        """
        Gets a :py:class:`NetworkEdge` from
        :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param net_cx: Network to extract edge from
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param edge_id: id of edge in `net_cx` network
        :type edge_id: int
        :return: Edge or `None` if not found
        :rtype: :py:class:`NetworkEdge`
        """
        edge_obj = net_cx.get_edge(edge_id)
        if edge_obj is None or edge_obj is (None, None):
            return None

        attributes = []
        e_attrs = net_cx.get_edge_attributes(edge_id)
        if e_attrs is not None:
            for e_attr in e_attrs:
                if 'd' in e_attr:
                    data_type = e_attr['d']
                else:
                    data_type = None
                e_annot = NetworkEdgeAttribute(name=e_attr['n'],
                                               value=e_attr['v'],
                                               data_type=data_type)
                attributes.append(e_annot)

        if 'i' in edge_obj:
            interaction = edge_obj['i']
        else:
            interaction = None
        s_node_obj = net_cx.get_node(edge_obj['s'])
        t_node_obj = net_cx.get_node(edge_obj['t'])
        n_edge = NetworkEdge(edge_id=edge_id, source_node_id=edge_obj['s'],
                             source_node_name=s_node_obj['n'],
                             target_node_id=edge_obj['t'],
                             target_node_name=t_node_obj['n'],
                             interaction=interaction,
                             attributes=attributes)
        return n_edge

    def get_all_edges_connected_to_node(self, net_cx=None, node_id=None):
        """
        Gets all edges connected to node
        :param net_cx:
        :param node_id:
        :return:
        """
        edge_list = []
        for edge_id, edge_obj in net_cx.get_edges():
            if edge_obj['s'] == node_id or edge_obj['t'] == node_id:
                net_edge = self.get_network_edge_from_network(net_cx=net_cx,
                                                              edge_id=edge_id)
                edge_list.append(net_edge)
        return edge_list

