""" Module to collapse PNVNode if they have the same sequence """
from __future__ import annotations
from typing import Dict, Set
from moPepGen import get_equivalent
from moPepGen.svgraph.PVGNode import PVGNode


class PVGCollapseNode(PVGNode):
    """ Node for collapsing """
    @classmethod
    def from_pvg_node(cls, node:PVGNode) -> PVGCollapseNode:
        """ Convert from a PVGNode """
        collapse_node = node.copy(in_nodes=False, out_nodes=True)
        collapse_node.__class__ = cls
        return collapse_node

    def __eq__(self, other:PVGCollapseNode):
        """ equal to """
        return self.seq.seq == other.seq.seq and \
            self.out_nodes == other.out_nodes and \
            self.cleavage == other.cleavage and \
            self.truncated == other.truncated and \
            self.reading_frame_index == other.reading_frame_index and \
            self.was_bridge == other.was_bridge

    def __ne__(self, other:PVGCollapseNode):
        """ not equal to """
        return not self == other

    def __hash__(self):
        """ hash """
        return hash((self.seq.seq, frozenset(self.out_nodes), self.cleavage,
            self.truncated, self.reading_frame_index, self.was_bridge))

class PVGNodeCollapser():
    """ Collapse PVGNode """
    def __init__(self, pool:Set[PVGNode]=None, mapper:Dict[PVGNode, PVGNode]=None):
        """ constructor """
        self.pool = pool or set()
        self.mapper = mapper or {}

    def collapse(self, node:PVGNode) -> PVGNode:
        """ Collapse the given node if it has the same in the pool """
        collapse_node = PVGCollapseNode.from_pvg_node(node)
        same_collapse_node = get_equivalent(self.pool, collapse_node)
        if same_collapse_node:
            same_node = self.mapper[same_collapse_node]
            if node.is_less_mutated(same_node):
                same_node.transfer_in_nodes_to(node)
                self.pool.remove(same_collapse_node)
                self.pool.add(collapse_node)
                self.mapper.pop(same_collapse_node, None)
                self.mapper[collapse_node] = node
                return same_node
            node.transfer_in_nodes_to(same_node)
            return node
        self.pool.add(collapse_node)
        self.mapper[collapse_node] = node
        return None