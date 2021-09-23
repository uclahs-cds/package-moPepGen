""" Module for node in the peptide graph """
from __future__ import annotations
import copy
from typing import List, Set
from moPepGen import aa, seqvar
from moPepGen.SeqFeature import FeatureLocation


class PVGNode():
    """ The PVGNode class is used in the svgraph.PeptideVariantGraph, that
    stores the sequence in a Node.

    Attributes:
        seq (aa.AminoAcidSeqRecord): The amino acid sequence.
        variants (List[seqvar.VariantRecordWithCoordinate]): The variant records
            carried in the sequence.
        frameshifts (Set[seqvar.VariantRecord]): Frameshifting variants.
        in_nodes (Set[PVGNode]): Inbound nodes.
        out_nodes (Set[PVGNode]): Outbound nodes
        cleavage (bool): Whether the start of the node is a cleavage site.
    """
    def __init__(self, seq:aa.AminoAcidSeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            in_nodes:Set[PVGNode]=None, out_nodes:Set[PVGNode]=None,
            frameshifts:Set[seqvar.VariantRecord]=None,
            cleavage:bool=False, truncated:bool=False, orf:List[int]=None):
        """ Construct a PVGNode object.

        Args:
            seq (aa.AminoAcidSeqRecordWithCoordinates): The amino acid sequence
                with the coordinates of that ORF (offset to the first M)
            variants (List[seqvar.VariantRecordWithCoordinate]): The variant records
                carried in the sequence.
            frameshifts (Set[seqvar.VariantRecord]): Frameshifting variants.
            in_nodes (Set[PVGNode]): Inbound nodes.
            ou_nodes (Set[PVGNode]): Outbound nodes
            cleavage (bool): Whether the start of the node is a cleavage site.
            truncated (bool): Whether the node is truncated. Useful when the
                sequence does not have a confirmed stop codon.
        """
        self.seq = seq
        self.variants = variants or []
        self.in_nodes = in_nodes or set()
        self.out_nodes = out_nodes or set()
        self.frameshifts = frameshifts or set ()
        self.cleavage = cleavage
        self.truncated = truncated
        self.orf = orf or [None, None]

    def __getitem__(self, index) -> PVGNode:
        """ get item """
        start, stop, _ = index.indices(len(self.seq))
        location = FeatureLocation(start=start, end=stop)
        seq = self.seq.__getitem__(index)
        variants = []
        frameshifts = copy.copy(self.frameshifts)

        for variant in self.variants:
            if variant.location.overlaps(location):
                variants.append(variant.shift(-start))
            elif variant.location.start >= stop and variant.varint in frameshifts:
                frameshifts.remove(variant.variant)
        return PVGNode(
            seq=seq,
            variants=variants,
            frameshifts=frameshifts,
            cleavage=self.cleavage,
            orf=self.orf
        )

    def add_out_edge(self, node:PVGNode) -> None:
        """ Add a outbound edge from this node.

        Args:
            node (PVGNode): The outboud node of the edge to add.
        """
        self.out_nodes.add(node)
        node.in_nodes.add(self)

    def remove_out_edge(self, node:PVGNode) -> None:
        """ Remove a outbound edge. It tries to remove the node, but won't
        raise any error if it doesn't exist.

        Args:
            node (PVGNode): The outbound node of the edge to remove.
        """
        try:
            self.out_nodes.remove(node)
        except KeyError:
            pass
        try:
            node.in_nodes.remove(self)
        except KeyError:
            pass

    def remove_out_edges(self) -> None:
        """ remove output nodes """
        for node in copy.copy(self.out_nodes):
            self.remove_out_edge(node)

    def find_reference_next(self) -> PVGNode:
        """ Find and return the next reference node. The next reference node
        is defined as the out node that has not variant, or not any variant
        that is not in this current node. """
        if not self.out_nodes:
            return None
        this_variants = [x.variant for x in self.variants]
        for node in self.out_nodes:
            if not node.variants:
                return node
            that_variants = [x.variant for x in node.variants]
            if not any(it not in this_variants for it in that_variants):
                return node
        raise ValueError('None of the out nodes is reference.')

    def find_reference_prev(self) -> PVGNode:
        """ Find and return the previous reference nocd. The previous reference
        node is defined as the inbond node that has not variant, or any variant
        that is not in this current node. """
        if not self.in_nodes:
            return None
        this_variants = [x.variant for x in self.variants]
        for node in self.in_nodes:
            if not node.variants:
                return node
            that_variants = [x.variant for x in node.variants]
            if not any(it not in this_variants for it in that_variants):
                return node
        raise ValueError('None of the in nodes is reference.')

    def split_node(self, index:int, cleavage:bool=False) -> PVGNode:
        """ Split the sequence at the given position, and create a new node
        as the outbound edge. Variants will also be adjusted. For example:

        KAPGCFP -> split at position 3 -> KAP-GCFP

        The node with GCFP is returned.

        Args:
            index (int): the position to split

        Returns:
            The new created node with the right part of the original sequence.
        """
        left_seq = self.seq[:index]
        right_seq = self.seq[index:]

        left_variants = []
        right_variants = []
        variants_to_pop = set()
        for variant in self.variants:
            if variant.location.start < index:
                left_variants.append(variant)
            if variant.location.end > index:
                right_variants.append(variant.shift(-index))
                if variant.location.start >= index:
                    variants_to_pop.add(variant.variant)
        self.seq = left_seq
        self.variants = left_variants

        new_node = PVGNode(seq=right_seq, variants=right_variants)
        new_node.frameshifts = copy.copy(self.frameshifts)
        new_node.orf = self.orf

        if cleavage:
            new_node.cleavage = True

        # need to remove the frameshift variant if it is no longer ther after
        # splitting the node.
        left_frameshifts = copy.copy(self.frameshifts)
        for variant in left_frameshifts:
            if variant in variants_to_pop:
                self.frameshifts.remove(variant)

        while self.out_nodes:
            node = self.out_nodes.pop()
            self.remove_out_edge(node)
            new_node.add_out_edge(node)
        self.add_out_edge(new_node)
        return new_node

    def truncate_right(self, i:int) -> PVGNode:
        """ Truncate the right i nucleotides off. """
        right_seq = self.seq[i:]
        left_variants = []
        right_variants = []
        right_frameshifts = copy.copy(self.frameshifts)

        for variant in self.variants:
            if variant.location.start < i:
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))
                if variant.location.start >= i and \
                        variant.variant.is_frameshifting():
                    self.frameshifts.remove(variant.variant)

        right_node = PVGNode(
            seq=right_seq,
            variants=right_variants,
            frameshifts=right_frameshifts,
        )

        self.seq = self.seq[:i]
        self.variants = left_variants

        return right_node

    def truncate_left(self, i:int) -> PVGNode:
        """ Truncate the left i nucleotides off. A new node with the left part
        of the sequences and variants associated is returned. The self node
        is updated with only the right part of the sequence and variants. """
        left_seq = self.seq[:i]
        left_variants = []
        right_variants = []
        left_frameshifts = copy.copy(self.frameshifts)

        for variant in self.variants:
            if variant.location.start < i:
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))
                if variant.location.start >= i and \
                        variant.variant.is_frameshifting():
                    left_frameshifts.remove(variant.variant)

        left_node = PVGNode(
            seq=left_seq,
            variants=left_variants,
            frameshifts=left_frameshifts
        )

        self.seq = self.seq[i:]
        self.variants = right_variants

        return left_node

    def find_start_index(self) -> int:
        """ Find the start amino acid position """
        return self.seq.seq.find('M')

    def copy(self) -> PVGNode:
        """ Create a copy of the node """
        return PVGNode(
            seq=self.seq,
            variants=copy.copy(self.variants),
            in_nodes=self.in_nodes,
            out_nodes=self.out_nodes,
            frameshifts=copy.copy(self.frameshifts),
            cleavage=self.cleavage,
            truncated=self.truncated,
            orf=self.orf
        )

    def get_nearest_next_ref_index(self) -> int:
        """ get the nearest reference index """
        cur = self
        while not cur.seq.locations and cur.out_nodes:
            cur = cur.find_reference_next()

        if not cur.seq.locations:
            return -1
        return int(cur.seq.locations[0].ref.start)

    def get_orf_start(self, i:int=0) -> int:
        """ get orf start position """
        for loc in self.seq.locations:
            if i < loc.query.end:
                return max(loc.query.start, i) * 3
        x = self.get_nearest_next_ref_index()
        if x == -1:
            return -1
        return x * 3
