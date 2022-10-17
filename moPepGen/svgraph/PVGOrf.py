""" PVGOrf """
from __future__ import annotations
from typing import Set, List
import copy
from moPepGen import circ, seqvar
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.svgraph.SubgraphTree import SubgraphTree


class PVGOrf():
    """ Helper class for an ORF and its corresponding start gain variants. """
    def __init__(self, orf:List[int,int]=None,
            start_gain:Set[seqvar.VariantRecord]=None, start_node:PVGNode=None):
        """ constructor """
        self.orf = orf or None
        self.start_gain = start_gain or set()
        self.start_node = start_node

    def copy(self) -> PVGOrf:
        """ copy """
        orf = self.orf
        start_gain = copy.copy(self.start_gain)
        return self.__class__(orf, start_gain, self.start_node)

    def __eq__(self, other:PVGOrf) -> bool:
        """ equal """
        return self.orf == other.orf and self.start_gain == other.start_gain

    def __ne__(self, other:PVGOrf) -> bool:
        """ not equal """
        return not self == other

    def __gt__(self, other:PVGOrf) -> bool:
        """ Greater than. We want it to be "greater" if the ORF is larger and
        less start gain variants or larger variant position. """
        if self.orf[0] > other.orf[0]:
            return True
        if self.orf[0] < other.orf[0]:
            return False
        if len(self.start_gain) < len(other.start_gain):
            return True
        if len(self.start_gain) > len(other.start_gain):
            return False
        for x, y in zip(self.start_gain, other.start_gain):
            if x > y:
                return True
            if x < y:
                return False
        return False

    def __ge__(self, other:PVGOrf) -> bool:
        """ larger then or equal to """
        return self > other or self == other

    def __lt__(self, other:PVGOrf) -> bool:
        """ less than """
        return not self >= other

    def __le__(self, other:PVGOrf) -> bool:
        """ less then or equal to """
        return not self > other

    def __hash__(self):
        """ hash """
        return hash((tuple(self.orf), frozenset(self.start_gain)))

    def is_valid_orf(self, node:PVGNode, subgraphs:SubgraphTree,
            circ_rna:circ.CircRNAModel) -> bool:
        """ """
        if not node.is_at_least_one_loop_downstream(self.start_node, subgraphs, circ_rna):
            return True
        start_gain = {x for x in self.start_gain if not x.is_circ_rna()}
        start_gain.update(x.variant for x in self.start_node.variants
            if not x.variant.is_circ_rna())
        variants = {x.variant for x in node.variants if not x.variant.is_circ_rna()}
        return not any(node.is_missing_variant(v) for v in start_gain) \
            and not any(x not in start_gain for x in variants)


