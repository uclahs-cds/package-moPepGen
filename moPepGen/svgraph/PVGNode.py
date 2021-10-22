""" Module for node in the peptide graph """
from __future__ import annotations
import copy
from functools import cmp_to_key
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
            reading_frame_index:int,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            in_nodes:Set[PVGNode]=None, out_nodes:Set[PVGNode]=None,
            frameshifts:Set[seqvar.VariantRecord]=None,
            cleavage:bool=False, truncated:bool=False, orf:List[int]=None,
            was_bridge:bool=False):
        """ Construct a PVGNode object.

        Args:
            seq (aa.AminoAcidSeqRecordWithCoordinates): The amino acid sequence
                with the coordinates reading frame
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
        self.frameshifts = frameshifts or set()
        self.cleavage = cleavage
        self.truncated = truncated
        self.orf = orf or [None, None]
        self.reading_frame_index = reading_frame_index
        self.was_bridge = was_bridge

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
            elif variant.location.start >= stop and variant.variant.is_frameshifting():
                frameshifts.discard(variant.variant)
        return PVGNode(
            seq=seq,
            variants=variants,
            frameshifts=frameshifts,
            cleavage=self.cleavage,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge
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
        self.out_nodes.discard(node)
        node.in_nodes.discard(self)

    def remove_out_edges(self) -> None:
        """ remove output nodes """
        for node in copy.copy(self.out_nodes):
            self.remove_out_edge(node)

    def is_bridge(self) -> None:
        """ Check if this is a bridge node to another reading frame """
        for node in self.out_nodes:
            if node.out_nodes and \
                    node.reading_frame_index != self.reading_frame_index:
                return True
        if self.was_bridge:
            return True
        return False

    def has_any_in_bridge(self) -> None:
        """ Check if it has any incoming node that is bridge """
        return any(node.is_bridge() for node in self.in_nodes)

    def get_variants_at(self, start:int, end:int=-1) -> seqvar.VariantRecord:
        """ Get the variant at position i """
        if end == -1:
            end = len(self.seq)
        if not -1 < start <= end <= len(self.seq):
            raise ValueError('start and end out of the range')
        variants = []
        location = FeatureLocation(start=start, end=end)
        for variant in self.variants:
            if variant.location.overlaps(location):
                variants.append(variant.variant)
        return variants

    def get_cleavage_gain_variants(self) -> List[seqvar.VariantRecord]:
        """ Get cleavage gain variants """
        cleavage_gain = []
        for variant in self.variants:
            if variant.location.end == len(self.seq):
                cleavage_gain.append(variant.variant)
        return cleavage_gain


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
        """ Find and return the previous reference node. The previous reference
        node is defined as the inbond node that has not variant, or any variant
        that is not in this current node. """
        if not self.in_nodes:
            return None
        if len(self.in_nodes) == 1:
            return list(self.in_nodes)[0]
        this_variants = [x.variant for x in self.variants]
        for node in self.in_nodes:
            if not node.variants:
                return node
            that_variants = [x.variant for x in node.variants]
            if not any(it not in this_variants for it in that_variants):
                return node
        raise ValueError('None of the in nodes is reference.')

    def find_least_varaint_prev(self) -> PVGNode:
        """ Find and return the inbonding node that has the least number of
        variants """
        if not self.in_nodes:
            return None
        if len(self.in_nodes) == 1:
            return list(self.in_nodes)[0]

        def sort_func(this:PVGNode, that:PVGNode):
            x = this.variants
            y = that.variants
            if len(x) < len(y):
                return -1
            if len(x) > len(y):
                return 1
            if len(x) == 0 and len(y) == 0:
                return -1
            for i,j in zip(x,y):
                if i.location < j.location:
                    return -1
                if i.location > j.location:
                    return 1
            return -1

        return sorted(self.in_nodes, key=cmp_to_key(sort_func))[0]

    def find_least_variant_next(self) -> PVGNode:
        """ Find the outbond node that has the least number of variants """
        if not self.out_nodes:
            return None
        if len(self.out_nodes) == 1:
            return list(self.out_nodes)[0]

        try:
            return self.find_reference_next()
        except ValueError:
            pass

        self_variants = {v.variant for v in self.variants}
        def sort_func(this:PVGNode, that:PVGNode):
            x = [v for v in this.variants if v.variant not in self_variants]
            y = [v for v in that.variants if v.variant not in self_variants]
            if len(x) < len(y):
                return -1
            if len(x) > len(y):
                return 1
            if len(x) == 0 and len(y) == 0:
                return -1
            for i,j in zip(x,y):
                if i.location < j.location:
                    return -1
                if i.location > j.location:
                    return 1
            return -1

        return sorted(self.out_nodes, key=cmp_to_key(sort_func))[0]

    def split_node(self, index:int, cleavage:bool=False) -> PVGNode:
        """ Split the sequence at the given position, and create a new node
        as the outbound edge. Variants will also be adjusted. For example:

        KAPGCFP -> split at position 3 -> KAP-GCFP

        The node with GCFP is returned.

        Args:
            index (int): the position to split
            cleavage (bool): If true, the cleavage attribute of the right node
                is update to be True.

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

        new_node = PVGNode(
            seq=right_seq,
            reading_frame_index=self.reading_frame_index,
            variants=right_variants,
            orf=self.orf,
            was_bridge=self.was_bridge
        )
        new_node.frameshifts = copy.copy(self.frameshifts)
        new_node.orf = self.orf

        if cleavage:
            new_node.cleavage = True

        # we only keep the last node to be was_bridge
        self.was_bridge = False

        # need to remove the frameshift variant if it is no longer ther after
        # splitting the node.
        for variant in variants_to_pop:
            self.frameshifts.discard(variant)

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
                    self.frameshifts.discard(variant.variant)

        right_node = PVGNode(
            seq=right_seq,
            variants=right_variants,
            frameshifts=right_frameshifts,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge
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
                    left_frameshifts.discard(variant.variant)

        left_node = PVGNode(
            seq=left_seq,
            variants=left_variants,
            frameshifts=left_frameshifts,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge
        )

        self.seq = self.seq[i:]
        self.variants = right_variants

        return left_node

    def append_left(self, other:PVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq

        frameshifts = copy.copy(other.frameshifts)
        for variant in self.variants:
            if variant.variant.is_frameshifting():
                frameshifts.add(variant.variant)
        self.frameshifts = frameshifts

        variants = copy.copy(other.variants)
        for variant in self.variants:
            variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants

    def append_right(self, other:PVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq

        for variant in other.variants:
            self.variants.append(variant.shift(len(self.seq.seq)))
            if variant.variant.is_frameshifting():
                self.frameshifts.add(variant)

        self.seq = new_seq

    def find_start_index(self) -> int:
        """ Find the start amino acid position """
        return self.seq.seq.find('M')

    def copy(self, in_nodes:bool=True, out_nodes:bool=True) -> PVGNode:
        """ Create a copy of the node """
        new_in_nodes = copy.copy(self.in_nodes) if in_nodes else set()
        new_out_nodes = copy.copy(self.out_nodes) if out_nodes else set()
        return PVGNode(
            seq=self.seq,
            variants=copy.copy(self.variants),
            in_nodes=new_in_nodes,
            out_nodes=new_out_nodes,
            frameshifts=copy.copy(self.frameshifts),
            cleavage=self.cleavage,
            truncated=self.truncated,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge
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
        """ get the ORF start position at the reference transcript sequence
        given M is found at the query postion i.
        """
        k = self.reading_frame_index
        if self.seq.locations:
            for loc in self.seq.locations:
                if i < loc.query.start:
                    return loc.ref.start * 3 + k
                if loc.query.start <= i < loc.query.end:
                    return (i - loc.query.start + loc.ref.start) * 3 + k

        out_node = self.find_least_variant_next()
        if str(out_node.seq.seq) == '*' and not out_node.out_nodes:
            return -1
        if out_node.seq.locations:
            return out_node.seq.locations[0].ref.start * 3 + k
        # raise ValueError('Can not find ORF')
        return -1
