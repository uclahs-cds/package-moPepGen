""" Module to collapse PNVNode if they have the same sequence """
from __future__ import annotations
import copy
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
        self_is_stop_altering = any(v.is_stop_altering for v in self.variants)
        other_is_stop_altering = any(v.is_stop_altering for v in other.variants)
        result = self.seq.seq == other.seq.seq \
            and not (self.seq.seq == '*' and self_is_stop_altering != other_is_stop_altering) \
            and self.out_nodes == other.out_nodes \
            and self.cleavage == other.cleavage \
            and self.truncated == other.truncated \
            and self.reading_frame_index == other.reading_frame_index \
            and self.was_bridge == other.was_bridge \
            and self.npop_collapsed == other.npop_collapsed == False \
            and {v.variant for v in self.variants if v.variant.is_alternative_splicing()} \
                == {v.variant for v in other.variants if v.variant.is_alternative_splicing()} \
            and {v.variant for v in self.variants if v.variant.is_indel()} \
                == {v.variant for v in other.variants if v.variant.is_indel()} \
            and self.subgraph_id == other.subgraph_id \
            and self.level == other.level \
            and (
                (
                    not any(v.variant.is_circ_rna() for v in self.variants)
                    and not any(v.variant.is_circ_rna() for v in other.variants)
                )
                or {
                    (v.variant.id, v.location.seqname)
                    for v in self.variants if not v.variant.is_circ_rna()
                } == {
                    (v.variant.id, v.location.seqname)
                    for v in other.variants if not v.variant.is_circ_rna()
                }
            )

        if result and hasattr(other, 'match'):
            other.match = self
        return result

    def __ne__(self, other:PVGCollapseNode):
        """ not equal to """
        return not self == other

    def __hash__(self):
        """ hash """
        alt_vars = {v.variant for v in self.variants if v.variant.is_alternative_splicing()}
        hash_values = (
            self.seq.seq, frozenset(self.out_nodes), self.cleavage,
            self.truncated, self.reading_frame_index, self.npop_collapsed,
            self.was_bridge, frozenset(alt_vars), self.subgraph_id, self.level
        )
        return hash(hash_values)

class PVGNodeCollapser():
    """ Collapse PVGNode """
    def __init__(self, pool:Set[PVGNode]=None, mapper:Dict[PVGNode, PVGNode]=None,
            is_circ_rna:bool=False):
        """ constructor """
        self.pool = pool or set()
        self.mapper = mapper or {}
        self.is_circ_rna = is_circ_rna

    @staticmethod
    def first_is_stop_altering(first:PVGNode, second:PVGNode) -> bool:
        """ Check whether the first node is stop altering while the second
        is not. """
        return any(x.is_stop_altering for x in first.variants) \
            and not any(y.is_stop_altering for y in second.variants)

    def should_keep_first(self, first:PVGNode, second:PVGNode) -> bool:
        """ Here we keep the node with least number of variants, unless one
        has stop altering mutation. """
        if self.first_is_stop_altering(first=first, second=second):
            return True
        if self.first_is_stop_altering(first=second, second=first):
            return False
        if first.is_less_mutated_than(second):
            return True
        return False

    def collapse(self, node:PVGNode) -> PVGNode:
        """ Collapse the given node if it has the same in the pool """
        collapse_node = PVGCollapseNode.from_pvg_node(node)
        same_collapse_node = get_equivalent(self.pool, collapse_node)
        original_upstreams = node.get_in_nodes()
        # If it already has the upstream indel map, the node must be already
        # collapsed at least once.
        upstream_indel_resolved = bool(node.upstream_indel_map)
        if same_collapse_node:
            same_node = self.mapper[same_collapse_node]

            if self.should_keep_first(node, same_node):
                same_node.transfer_in_nodes_to(node)
                self.pool.remove(same_collapse_node)
                self.pool.add(collapse_node)
                self.mapper.pop(same_collapse_node, None)
                self.mapper[collapse_node] = node
                node_to_keep = node
                node_to_discard = same_node
            else:
                node.transfer_in_nodes_to(same_node)
                node_to_keep = same_node
                node_to_discard = node

            for out_node in node_to_discard.get_out_nodes():
                if node_to_discard in out_node.upstream_indel_map:
                    indels = out_node.upstream_indel_map.pop(node_to_discard)
                    out_node.upstream_indel_map[node_to_keep] = indels

            indel_map = {k:copy.copy(v) for k,v in node_to_discard.upstream_indel_map.items()}
            node_to_keep.upstream_indel_map.update(indel_map)
        else:
            self.pool.add(collapse_node)
            self.mapper[collapse_node] = node
            node_to_keep = node
            node_to_discard = None

        if not upstream_indel_resolved and node.has_any_indel():
            indels = [x.variant for x in node.variants
                if x.variant.is_indel() and not x.downstream_cleavage_altering]
            for upstream in original_upstreams:
                node_to_keep.upstream_indel_map[upstream] = copy.copy(indels)

        return node_to_discard

class PVGPopCollapseNode(PVGNode):
    """ Node for collapsing """
    @classmethod
    def from_pvg_node(cls, node:PVGNode) -> PVGPopCollapseNode:
        """ Convert from a PVGNode """
        pop_collapse_node = node.copy(in_nodes=False, out_nodes=True)
        pop_collapse_node.__class__ = cls
        return pop_collapse_node

    def __eq__(self, other:PVGPopCollapseNode):
        """ equal to """
        self_is_stop_altering = any(v.is_stop_altering for v in self.variants)
        other_is_stop_altering = any(v.is_stop_altering for v in other.variants)

        result = self.seq.seq == other.seq.seq \
            and not (self.seq.seq == '*' and self_is_stop_altering != other_is_stop_altering) \
            and self.out_nodes == other.out_nodes \
            and self.cleavage == other.cleavage \
            and self.truncated == other.truncated \
            and self.reading_frame_index == other.reading_frame_index \
            and self.was_bridge == other.was_bridge \
            and self.npop_collapsed == other.npop_collapsed \
            and self.cpop_collapsed == other.cpop_collapsed \
            and self.subgraph_id == other.subgraph_id \
            and {v.variant for v in self.variants if v.variant.is_alternative_splicing()}\
                == {v.variant for v in other.variants if v.variant.is_alternative_splicing()} \
            and {v.variant for v in self.variants if v.variant.is_indel()} \
                == {v.variant for v in other.variants if v.variant.is_indel()} \
            and self.subgraph_id == other.subgraph_id \
            and self.level == other.level \
            and (
                (
                    not any(v.variant.is_circ_rna() for v in self.variants)
                    and not any(v.variant.is_circ_rna() for v in other.variants)
                )
                or {
                    (v.variant.id, v.location.seqname)
                    for v in self.variants if not v.variant.is_circ_rna()
                } == {
                    (v.variant.id, v.location.seqname)
                    for v in other.variants if not v.variant.is_circ_rna()
                }
            )

        if result and hasattr(other, 'match'):
            other.match = self
        return result

    def __ne__(self, other:PVGPopCollapseNode):
        """ not equal to """
        return not self == other

    def __hash__(self):
        """ hash """
        alt_vars = {v.variant for v in self.variants if v.variant.is_alternative_splicing()}
        hash_items = (
            self.seq.seq, frozenset(self.out_nodes), self.cleavage,
            self.truncated, self.reading_frame_index, self.was_bridge,
            self.npop_collapsed, self.cpop_collapsed, self.subgraph_id,
            frozenset(alt_vars), self.subgraph_id, self.level
        )
        return hash(hash_items)

class PVGNodePopCollapser():
    """ Pop collapse PVGNodes. Similar to PVGNodeCollapser, but only nodes
    with smallest non-zero variants are kept unless all nodes are reference. """
    def __init__(self, pool:Set[PVGPopCollapseNode]=None,
            mapper:Dict[PVGPopCollapseNode, PVGNode]=None):
        """ constructor """
        self.pool = pool or set()
        self.mapper = mapper or {}

    def collapse(self, node:PVGNode) -> PVGNode:
        """ Collapse the given node if it has the same in the pool """
        collapse_node = PVGPopCollapseNode.from_pvg_node(node)
        same_collapse_node = get_equivalent(self.pool, collapse_node)
        if same_collapse_node:
            same_node = self.mapper[same_collapse_node]
            if node.is_less_mutated_than(same_node):
                same_node.transfer_in_nodes_to(node)
                self.pool.remove(same_collapse_node)
                self.pool.add(collapse_node)
                self.mapper.pop(same_collapse_node, None)
                self.mapper[collapse_node] = node
                node_to_keep = node
                node_to_discard = same_node
            else:
                node.transfer_in_nodes_to(same_node)
                node_to_keep = same_node
                node_to_discard = node

            for out_node in node_to_discard.get_out_nodes():
                if node_to_discard in out_node.upstream_indel_map:
                    indels = out_node.upstream_indel_map.pop(node_to_discard)
                    out_node.upstream_indel_map[node_to_keep] = indels

        else:
            self.pool.add(collapse_node)
            self.mapper[collapse_node] = node
            node_to_keep = node
            node_to_discard = None

        if node.has_any_indel():
            indels = [x.variant for x in node.variants
                if x.variant.is_indel() and not x.downstream_cleavage_altering]
            for upstream in node.in_nodes:
                node_to_keep.upstream_indel_map[upstream] = indels

        return node_to_discard
