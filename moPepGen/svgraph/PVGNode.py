""" Module for node in the peptide graph """
from __future__ import annotations
import copy
from typing import List, Set
from moPepGen import aa, seqvar


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
    def __init__(self, seq:aa.AminoAcidSeqRecord,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            in_nodes:Set[PVGNode]=None, out_nodes:Set[PVGNode]=None,
            frameshifts:Set[seqvar.VariantRecord]=None,
            cleavage:bool=False, truncated:bool=False, orf:List[int]=None):
        """ Construct a PVGNode object.

        Args:
            seq (aa.AminoAcidSeqRecord): The amino acid sequence.
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
        self.variants = variants
        self.in_nodes = set() if in_nodes is None else in_nodes
        self.out_nodes = set() if out_nodes is None else out_nodes
        self.frameshifts = set() if frameshifts is None else frameshifts
        self.cleavage = cleavage
        self.truncated = truncated
        self.orf = orf if orf is not None else [None, None]

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
        node is defined as the inbond node that has not variant, or all variant
        from the current node can be found in that node. """
        if not self.in_nodes:
            return None
        this_variants = [x.variant for x in self.variants]
        for node in self.in_nodes:
            if not node.variants:
                return node
            that_variants = [x.variant for x in node.variants]
            if not any(it not in that_variants for it in this_variants):
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
