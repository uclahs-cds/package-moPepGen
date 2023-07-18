""" Module to collapse TNGNode if they have the same sequence """
from __future__ import annotations
from typing import Dict, Set
from moPepGen import get_equivalent
from moPepGen.svgraph.TVGEdge import TVGEdge
from moPepGen.svgraph.TVGNode import TVGNode


class TVGCollapseNode(TVGNode):
    """ Node for collapsing """
    @classmethod
    def from_tvg_node(cls, node:TVGNode) -> TVGCollapseNode:
        """ Convert from a PVGNode """
        collapse_node = node.copy()
        collapse_node.__class__ = cls
        for edge in node.in_edges:
            edge_copy = TVGEdge(edge.in_node, collapse_node, edge.type)
            collapse_node.in_edges.add(edge_copy)
        for edge in node.out_edges:
            edge_copy = TVGEdge(collapse_node, edge.out_node, edge.type)
            collapse_node.out_edges.add(edge_copy)
        return collapse_node

    def __eq__(self, other:TVGCollapseNode):
        """ equal to """
        result = self.seq.seq == other.seq.seq \
            and set(self.get_in_nodes()) == set(other.get_in_nodes()) \
            and set(self.get_out_nodes()) == set(other.get_out_nodes()) \
            and self.reading_frame_index == other.reading_frame_index \
            and self.subgraph_id == other.subgraph_id \
            and self.was_bridge == other.was_bridge

        if result and hasattr(other, 'match'):
            other.match = self
        return result

    def __ne__(self, other:TVGCollapseNode):
        """ not equal to """
        return not self == other

    def __hash__(self):
        """ hash """
        hash_values = (
            self.seq.seq, frozenset(self.get_out_nodes()),
            frozenset(self.get_in_nodes()), self.reading_frame_index,
            self.subgraph_id, self.was_bridge
        )
        return hash(hash_values)

class TVGNodeCollapser():
    """ Collapse PVGNode """
    def __init__(self, pool:Set[TVGNode]=None, mapper:Dict[TVGNode, TVGNode]=None):
        """ constructor """
        self.pool = pool or set()
        self.mapper = mapper or {}

    @staticmethod
    def first_is_stop_altering(first:TVGNode, second:TVGNode) -> bool:
        """ Check whether the first node is stop altering while the second
        is not. """
        return any(x.is_stop_altering for x in first.variants) \
            and not any(y.is_stop_altering for y in second.variants)

    def should_keep_first(self, first:TVGNode, second:TVGNode) -> bool:
        """ Here we keep the node with least number of variants, unless one
        has stop altering mutation. """
        if self.first_is_stop_altering(first=first, second=second):
            return True
        if self.first_is_stop_altering(first=second, second=first):
            return False
        if first.is_less_mutated_than(second):
            return True
        return False

    def collapse(self, node:TVGNode) -> TVGNode:
        """ Collapse the given node if it has the same in the pool """
        collapse_node = TVGCollapseNode.from_tvg_node(node)
        same_collapse_node = get_equivalent(self.pool, collapse_node)
        if same_collapse_node:
            same_node = self.mapper[same_collapse_node]

            if self.should_keep_first(node, same_node):
                self.pool.remove(same_collapse_node)
                self.pool.add(collapse_node)
                self.mapper.pop(same_collapse_node, None)
                self.mapper[collapse_node] = node
                node_to_discard = same_node
            else:
                node_to_discard = node
        else:
            self.pool.add(collapse_node)
            self.mapper[collapse_node] = node
            node_to_discard = None

        return node_to_discard
