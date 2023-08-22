""" PVGOrf """
from __future__ import annotations
from typing import Set, List
import copy
from moPepGen import circ, seqvar
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.svgraph.SubgraphTree import SubgraphTree


class PVGOrf():
    """ Helper class for an ORF and its corresponding start gain variants. """
    def __init__(self, orf:List[int,int]=None,
            start_gain:Set[seqvar.VariantRecord]=None, start_node:PVGNode=None,
            subgraph_id:str=None, node_offset:int=None,
            cleavage_gain:Set[seqvar.VariantRecord]=None):
        """ constructor """
        self.orf = orf or None
        self.start_gain = start_gain or set()
        self.start_node = start_node
        self.subgraph_id = subgraph_id
        self.node_offset = node_offset
        self._orf_node_locs = None
        self.cleavage_gain = cleavage_gain or set()

    def copy(self) -> PVGOrf:
        """ copy """
        orf = self.orf
        start_gain = copy.copy(self.start_gain)
        cleavage_gain = copy.copy(self.cleavage_gain)
        return self.__class__(
            orf=orf, start_gain=start_gain, start_node=self.start_node,
            subgraph_id=self.subgraph_id, node_offset=self.node_offset,
            cleavage_gain=cleavage_gain
        )

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

    def get_orf_node_locs(self) -> List[FeatureLocation]:
        """ Get locations of the orf node. """
        if self._orf_node_locs:
            return self._orf_node_locs
        locs = {
            loc.get_ref_dna_location() for loc in self.start_node.seq.locations
        }
        locs.update([
            v.variant.location for v in self.start_node.variants
            if not v.variant.is_circ_rna()
        ])
        self._orf_node_locs = locs
        return sorted(locs)

    def location_is_at_least_one_loop_downstream(self, i:int, subgraph_id:str,
            subgraphs:SubgraphTree, circ_rna:circ.CircRNAModel) -> bool:
        """ Whether the genetic or transcriptional location i with the given
        subgraph ID is at least one loop downstream to the start of the node
        where the ORF start is located.

        Args:
            - `i` (int): The genetic or transcriptional (depend on the coordinate
              used by the graph) downstream to the ORF to check.
            - `subgraph_id` (str): The subgraph ID of the node where the
              location `i` is located.
            - `subgraphs` (SubgraphTree): The subgraph tree object which contains
              the relationship of each subgraph. Should be passed from the graph.
            - `circ_rna` (CircRNAModel): The circRNA record. Its fragments are
              used to tell which position appears earlier in the sequence of
              the circRNA.
        """
        x_level = subgraphs[subgraph_id].level

        if self.start_node.seq.locations \
                and 0 in self.start_node.seq.locations[0].query:
            loc = self.start_node.seq.locations[0]
            y_level = subgraphs[loc.ref.seqname].level
            j = loc.get_ref_dna_start()
        else:
            for v in self.start_node.variants:
                if v.variant.is_circ_rna():
                    continue
                if 0 in v.location:
                    j = v.variant.location.start
                    y_level = subgraphs[v.location.seqname].level
                    break
            else:
                raise ValueError('Failed to find the position for start of ORF.')

        if x_level - y_level > 1:
            return True

        if x_level == y_level:
            return False

        if x_level - y_level == 1:
            for fragment in circ_rna.fragments:
                frag = fragment.location
                if i in frag:
                    if j in frag:
                        return i >= j
                    if j + 1 in frag:
                        return i >= j + 1
                if i in frag:
                    return False
                if j in frag or j + 1 in frag:
                    return True
        raise ValueError('Locations not found from the fragments of the circRNA.')

    def node_is_at_least_one_loop_downstream(self, node:PVGNode, subgraphs:SubgraphTree,
            circ_rna:circ.CircRNAModel):
        """ Checks whether a given node is at least one loop downstream to the
        ORF start site. """
        if node.seq.locations:
            loc = node.seq.locations[-1]
            i = loc.get_ref_dna_end()
            subgraph_id = node.seq.locations[-1].ref.seqname
        else:
            for v in reversed(node.variants):
                if not v.variant.is_circ_rna():
                    i = v.variant.location.end - 1
                    subgraph_id = v.location.seqname
                    break
            else:
                raise ValueError("Failed to find a non circRNA variant.")

        return self.location_is_at_least_one_loop_downstream(i, subgraph_id,
            subgraphs, circ_rna)

    def is_valid_orf(self, node:PVGNode, subgraphs:SubgraphTree,
            circ_rna:circ.CircRNAModel,
            upstream_variants:Set[seqvar.VariantRecord]=None) -> bool:
        """ Checks if it is a valid orf of a downstream node. """
        if node.is_hybrid_node(subgraphs):
            return False
        upstream_variants = upstream_variants or set()

        if not self.node_is_at_least_one_loop_downstream(node, subgraphs, circ_rna):
            return not self.has_any_incompatible_variant(node)

        start_gain = {x for x in self.start_gain if not x.is_circ_rna()}
        start_gain.update(upstream_variants)
        start_gain.update(x.variant for x in self.start_node.variants
            if not x.variant.is_circ_rna() and x.not_cleavage_altering())
        start_gain = {v for v in start_gain if not v.is_circ_rna()}

        variants = set()
        for v in node.variants:
            if v.variant.is_circ_rna():
                continue
            i = v.variant.location.end - 1
            subgraph_id = v.location.seqname
            if self.location_is_at_least_one_loop_downstream(i, subgraph_id,
                    subgraphs, circ_rna):
                variants.add(v.variant)
            else:
                start_gain.add(v.variant)

        return not any(node.is_missing_variant(v) for v in start_gain) \
            and not any(x not in start_gain for x in variants)

    def is_valid_orf_to_misc_nodes(self, nodes:List[PVGNode], subgraphs:SubgraphTree,
            circ_rna:circ.CircRNAModel) -> bool:
        """ Checks if the ORF is valid for a given series of miscleaved peptide
        nodes. """
        upstream_variants = set()
        for node in nodes:
            if not self.is_valid_orf(node, subgraphs, circ_rna, upstream_variants):
                return False
            upstream_variants.update(x.variant for x in node.variants)

        downstream_valid = False
        for out_node in nodes[-1].get_out_nodes():
            if (out_node.seq.seq == '*' and not out_node.get_out_nodes()) \
                    or self.is_valid_orf(out_node, subgraphs, circ_rna, upstream_variants):
                downstream_valid = True

        return downstream_valid

    def has_any_incompatible_variant(self, node:PVGNode) -> bool:
        """ Checks if the orf is missing any frameshift variants before the
        orf start node. """
        for v in node.variants:
            if v.variant.is_circ_rna():
                continue

            if not self.is_compatible_with_variant(v):
                return True
        return False

    def is_compatible_with_variant(self, variant:seqvar.VariantRecordWithCoordinate) -> bool:
        """ Checks if the ORF is compatible with a given variant. A variant is
        incompatible if it overlaps with any of the locations of the ORF start
        node but not carried by the start node or is not a start gain variant
        already. """
        orf_variants = {v.variant.id for v in self.start_node.variants}
        orf_variants.update({v.id for v in self.start_gain})

        orf_node_locs = self.get_orf_node_locs()

        orf_subgraph = self.start_node.get_first_subgraph_id()

        return variant.location.seqname == orf_subgraph \
            or not any(variant.variant.location.overlaps(loc) for loc in orf_node_locs) \
            or variant.variant.id in orf_variants
