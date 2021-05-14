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
        self._equalityfailreason = None

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

    def __eq__(self, other):
        """
        Compares self with **other** for equality
        :param other:
        :type other: :py:class:`Attribute`
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

    def __str__(self):
        """
        Returns str representation of object
        :return:
        """
        return 'name=' + str(self._name) + ', values=' +\
               str(self._value) + ', type=' +\
               str(self._data_type)

    def get_equality_fail_reason(self):
        """

        :return:
        """
        return self._equalityfailreason

    def add_attribute_to_edge(self, net_cx=None, edge_id=None):
        """
        Adds this edge attribute

        :param net_cx: Network to add edge attribute to
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param edge_id: Id of edge on **net_cx** network
        :type edge_id: int
        :return:
        """
        net_cx.set_edge_attribute(edge_id, self._name,
                                  self._value,
                                  type=self._data_type)

    def add_attribute_to_node(self, net_cx=None, node_id=None):
        """
        Adds this node attribute

        :param net_cx: Network to add node attribute to
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param node_id: Id of node on **net_cx** network
        :type node_id: int
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
    Object that represents an edge in a Network

    """
    def __init__(self, edge_id=None, source_node_id=None,
                 source_node_name=None,
                 target_node_id=None,
                 target_node_name=None, interaction=None,
                 attributes=None):
        """
        Constructor

        :param edge_id: Id of edge
        :type edge_id: int
        :param source_node_id: Id of source node
        :type source_node_id: int
        :param source_node_name: Name of source node
        :type source_node_name: str
        :param target_node_id: Id of target node
        :type target_node_id: int
        :param target_node_name: Name of target node
        :type target_node_name: str
        :param interaction: Interaction for edge
        :type interaction: str
        :param attributes: Edge :py:class:`Attribute`
        :type attributes: list
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

        :return: Id of edge
        :rtype: int
        """
        return self._edge_id

    def set_id(self, edge_id):
        """
        Sets edge id

        :param edge_id: Id of edge
        :type edge_id: int
        """
        self._edge_id = edge_id

    def get_source_node_id(self):
        """
        Gets source node id

        :return: Id of source node
        :rtype: int
        """
        return self._source_node_id

    def set_source_node_id(self, source_node_id):
        """
        Sets source node id

        :param source_node_id: Id of source node
        :type source_node_id: int
        """
        self._source_node_id = source_node_id

    def get_target_node_id(self):
        """
        Gets target node id

        :return: Id of target node id
        :rtype: int
        """
        return self._target_node_id

    def set_target_node_id(self, target_node_id):
        """
        Sets target node id

        :param target_node_id: Id of target node
        :type target_node_id: int
        """
        self._target_node_id = target_node_id

    def get_source_node_name(self):
        """
        Gets source node name

        :return: Name of source node
        :rtype: str
        """
        return self._source_node_name

    def set_source_node_name(self, node_name):
        """
        Sets source node name

        :param node_name: Name of source node
        :type node_name: str
        """
        self._source_node_name = node_name

    def get_target_node_name(self):
        """
        Gets target node name

        :return: Target node name
        :rtype: str
        """
        return self._target_node_name

    def set_target_node_name(self, node_name):
        """
        Sets target node name

        :param node_name: Target node name
        :type node_name: str
        """
        self._target_node_name = node_name

    def get_interaction(self):
        """
        Gets interaction

        :return: Interaction
        :rtype: str
        """
        return self._interaction

    def set_interaction(self, interaction):
        """
        Sets interaction

        :param interaction: Interaction
        :type interaction: str
        """
        self._interaction = interaction

    def get_attributes(self):
        """
        Gets edge annotations

        :return: :py:class:`Attribute` as a list
        :rtype: list
        """
        return self._attributes

    def set_attributes(self, attributes):
        """
        Sets edge attributes, replacing any existing
        attributes

        :param attributes: list of :py:class:`Attribute`
        :type attributes: list
        """
        self._attributes = attributes

    def add_edge_to_network(self, net_cx=None, source_node_id=None,
                            target_node_id=None):
        """
        Adds this edge and its attributes to network **net_cx** ignoring
        the source and target node values in this object and using
        the ones passed in.

        :param net_cx: Network to modify
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param source_node_id: Id of source node to use for edge. To use
                               source node id of this object leave as
                               default of `None`
        :type source_node_id: int
        :param target_node_id: Id of target node to use for edge. To use
                               target node id of this object leave as
                               default of `None`
        :return: Id of new edge
        :rtype: int
        """
        if source_node_id is None:
            src_node = self._source_node_id
        else:
            src_node = source_node_id
        if target_node_id is None:
            target_node = self._target_node_id
        else:
            target_node = target_node_id
        new_edge_id = net_cx.create_edge(edge_source=src_node,
                                         edge_target=target_node,
                                         edge_interaction=self._interaction)
        if self._attributes is not None:
            for annot in self._attributes:
                annot.add_attribute_to_edge(net_cx=net_cx,
                                            edge_id=new_edge_id)
        return new_edge_id

    def remove_edge_from_network(self, net_cx=None):
        """
        Removes this edge from **net_cx** network passed in. Any
        attributes connected with this edge are also removed.

        :param net_cx: Network to modify
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
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
        Gets string representation of edge in this format

        .. code-block:: python

            po=<EDGE ID>, s=<SRC NODE ID> (<SRC NODE NAME>),\
             t=<TARGET NODE ID> (<TARGET NODE NAME>),\
             i=<INTERACTION>

        :return: String representation of edge ``None`` values will
                 be printed as ``None``
        :rtype: str
        """
        return 'po=' + str(self._edge_id) + ', s=' +\
               str(self._source_node_id) + ' (' +\
               str(self._source_node_name) + '), t=' +\
               str(self._target_node_id) + ' (' +\
               str(self._target_node_name) + '), i=' +\
               str(self._interaction)

    def merge_edge(self, other_edge):
        """
        Merges edge with this one. source, target, interaction is ignored
        and attributes are combined with duplicates removed. In case of
        attributes with matching names and different values those values
        are switched to a list of matching type and if types differ
        then list_of_string is used

        :param other_edge:
        :return:
        """
        if other_edge is None:
            return

        if other_edge.get_attributes() is None:
            return
        if self._attributes is None:
            self._attributes = []
        self._attributes.extend(other_edge.get_attributes())

        self.remove_identical_attributes()
        self.merge_duplicate_attributes()

    def merge_duplicate_attributes(self):
        """

        :return:
        """
        if self._attributes is None:
            return

        # build dict of lists with matching attribute names
        attrs_dict = {}
        for me_attr in self._attributes:
            if me_attr.get_name() not in attrs_dict:
                attrs_dict[me_attr.get_name()] = []
            attrs_dict[me_attr.get_name()].append(me_attr)

        new_attrs = []
        for key in attrs_dict.keys():
            # simple case only 1 attribute with name
            if len(attrs_dict[key]) == 1:
                new_attrs.append(attrs_dict[key][0])
                continue
            new_attrs.append(self._merge_attributes(attrs_dict[key]))
        self._attributes = new_attrs

    def _get_data_type_for_attributes(self, attrs_list):
        """
        Given a list of :py:class:`Attribute` objects examine
        the data types and return a data type based on these
        rules:


        If data types are same <type> or list_of_<type>
        just return list_of_<type>

        If data types are different but one is list_of_<type> and
        other is <type> then return list_of_<type>

        All other cases return list_of_string

        :param attrs_list:
        :return:
        """
        d_type_set = set()
        # build set of data types
        for attr in attrs_list:
            d_type_set.add(attr.get_data_type())

        # easier case only 1 data type
        if len(d_type_set) == 1:
            d_type = d_type_set.pop()

            # if string or None list_of_string
            if d_type is None or d_type == 'string':
                return 'list_of_string'

            # if <type> just return list_of_<type>
            if not d_type.startswith('list_of'):
                return 'list_of_' + d_type
            # already was set to list, leave it
            return d_type

        # harder case 2 data types
        if len(d_type_set) == 2:
            d_type_one = d_type_set.pop()
            d_type_two = d_type_set.pop()

            # check if data types are <type> and list_of_<type>
            # if so return list_of_<type>
            if d_type_one is not None and d_type_one.startswith('list_of_'):
                if d_type_two is not None and d_type_one.endswith(d_type_two):
                    return d_type_one
            if d_type_two is not None and d_type_two.startswith('list_of_'):
                if d_type_one is not None and d_type_two.endswith(d_type_one):
                    return d_type_two

        # failure case if none of the above worked
        # just go with list_of_string
        return 'list_of_string'

    def _merge_attributes(self, attrs_list):
        """
        Takes list of :py:class:`Attribute` objects and merges
        them following these rules.

        If data types are the same they are added to a set and that
        set is converted to a list and set as value for this attribute

        If data types differ, then all values are added as str and
        put in a set. This set is converted to a list and set
        as value for this attribute. The data type is then set
        to 'list_of_string'

        .. note::

            Attribute name is taken from 1st :py:class:`Attribute` in list

        :param attrs_list: list with more than 1 :py:class:`Attribute` objects
        :type attrs_list: list
        :return: Merged edge attribute
        :rtype: :py:class:`Attribute`
        """
        d_type = self._get_data_type_for_attributes(attrs_list)

        name = attrs_list[0].get_name()
        attr_value_set = set()
        for attr in attrs_list:
            if isinstance(attr.get_value(), list):
                for el in attr.get_value():
                    if d_type == 'list_of_string':
                        attr_value_set.add(str(el))
                    else:
                        attr_value_set.add(el)
            else:
                if d_type == 'list_of_string':
                    attr_value_set.add(str(attr.get_value()))
                else:
                    attr_value_set.add(attr.get_value())
        return Attribute(name=name,
                         value=list(attr_value_set),
                         data_type=d_type)

    def remove_identical_attributes(self):
        """
        Examines attributes contained within and
        keeps only one of duplicates.

        .. note::

            This method recreates the attributes list within

        :return: list of :py:class:`Attribute` objects that were
                 duplicates and removed
        :rtype list:
        """
        merges = []
        if self._attributes is None:
            return merges

        if len(self._attributes) <= 1:
            return []

        new_attrs = []
        me_attr = self._attributes.pop()
        while me_attr is not None:
            found_match = False
            for other_attr in self._attributes:
                if me_attr == other_attr:
                    merges.append(me_attr)
                    found_match = True
                    break
            if found_match is False:
                new_attrs.append(me_attr)
            if len(self._attributes) == 0:
                break
            me_attr = self._attributes.pop()

        self._attributes = new_attrs
        return merges


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
                n_annot = Attribute(name=n_attr['n'],
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
                e_annot = Attribute(name=e_attr['n'],
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

    def get_all_edges_connected_to_node(self, net_cx=None, node_id=None,
                                        interaction=None):
        """
        Gets all edges connected to node
        :param net_cx: Network to get edges from
        :type net_cx: :py:class:`~ndex2.nice_cx_network.NiceCXNetwork`
        :param node_id: Id of node
        :type node_id: int
        :param interaction: If set only include edges whose interaction matches
                            this value. Leaving as ``None`` will include all edges.
        :type interaction: str
        :return: List of :py:class:`NetworkEdge` objects connected to to **node_id**
        :rtype: list
        """
        edge_list = []
        for edge_id, edge_obj in net_cx.get_edges():
            if edge_obj['s'] == node_id or edge_obj['t'] == node_id:
                if interaction is None or edge_obj['i'] == interaction:
                    net_edge = self.get_network_edge_from_network(net_cx=net_cx,
                                                                  edge_id=edge_id)
                    edge_list.append(net_edge)
        return edge_list




