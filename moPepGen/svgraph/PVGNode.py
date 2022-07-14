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
        reading_frame_index (int): Reading frame index that the node belongs to.
            Must be 0, 1, or 2.
        variants (List[seqvar.VariantRecordWithCoordinate]): The variant records
            carried in the sequence.
        in_nodes (Set[PVGNode]): Inbound nodes.
        out_nodes (Set[PVGNode]): Outbound nodes
        cleavage (bool): Whether the start of the node is a cleavage site.
        truncated (bool): Whether the node is truncated. Useful when the
            sequence does not have a confirmed stop codon.
        orf (List[int]): List of two integers indicating the start and stop
            of the orf.
        was_bridge (bool): When true indicating that the node used to be a
            bridge node between reading frames, but is no longer a bridge node
            because it was split.
        pre_cleaved (bool): Indicating that the node is already cleaved.
        level (bool): Subgraph level that the node belongs to.
        npop_collapsed (bool): Whether the n-terminus is pop collapsed.
        cpop_collapsed (bool): Whether the c-terminus is pop collapsed.

    """
    def __init__(self, seq:aa.AminoAcidSeqRecordWithCoordinates,
            reading_frame_index:int, subgraph_id:str,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            in_nodes:Set[PVGNode]=None, out_nodes:Set[PVGNode]=None,
            cleavage:bool=False, truncated:bool=False, orf:List[int]=None,
            was_bridge:bool=False, pre_cleaved:bool=False, level:int=0,
            npop_collapsed:bool=False, cpop_collapsed:bool=False):
        """ Construct a PVGNode object. """
        self.seq = seq
        self.variants = variants or []
        self.in_nodes = in_nodes or set()
        self.out_nodes = out_nodes or set()
        self.cleavage = cleavage
        self.truncated = truncated
        self.orf = orf or [None, None]
        self.reading_frame_index = reading_frame_index
        self.was_bridge = was_bridge
        self.pre_cleave = pre_cleaved
        self.subgraph_id = subgraph_id
        self.level = level
        self.npop_collapsed = npop_collapsed
        self.cpop_collapsed = cpop_collapsed

    def __getitem__(self, index) -> PVGNode:
        """ get item """
        start, stop, _ = index.indices(len(self.seq))
        location = FeatureLocation(start=start, end=stop)
        seq = self.seq.__getitem__(index)
        variants = []

        for variant in self.variants:
            if variant.location.overlaps(location):
                variants.append(variant.shift(-start))

        npop_collapsed = self.npop_collapsed and start == 0
        cpop_collapsed = self.cpop_collapsed and stop == len(seq) or stop == -1

        return PVGNode(
            seq=seq,
            variants=variants,
            cleavage=self.cleavage,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            npop_collapsed=npop_collapsed,
            cpop_collapsed=cpop_collapsed
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

    def get_out_nodes(self) -> List[PVGNode]:
        """ Get outgoing nodes as a list. """
        return list(self.out_nodes)

    def get_in_nodes(self) -> List[PVGNode]:
        """ Get incoming nodes as a list """
        return list(self.in_nodes)

    def remove_out_edges(self) -> None:
        """ remove output nodes """
        for node in copy.copy(self.out_nodes):
            self.remove_out_edge(node)

    def is_inbond_of(self, node:PVGNode) -> bool:
        """ Checks if self is the inbond node of a given node. """
        return node in self.out_nodes

    def is_orphan(self) -> bool:
        """ Checks if the node is orphan (no inbond or outbond node) """
        return not self.in_nodes and not self.out_nodes

    def is_bridge(self) -> None:
        """ Check if this is a bridge node to another reading frame """
        for node in self.out_nodes:
            if not (not node.out_nodes and node.seq.seq == '*') and \
                    node.reading_frame_index != self.reading_frame_index:
                return True
        if self.was_bridge:
            return True
        return False

    def is_subgraph_end(self) -> bool:
        """ Check if this is the end node of a subgraph """
        return all(x.subgraph_id == self.subgraph_id for x in self.in_nodes) and \
            any(x.subgraph_id != self.subgraph_id for x in self.out_nodes)

    def is_subgraph_start(self) -> bool:
        """ Check if this is the start node of a subgraph """
        return any(x.subgraph_id != self.subgraph_id for x in self.in_nodes) and \
            all(x.subgraph_id == self.subgraph_id for x in self.out_nodes)

    def has_any_in_bridge(self) -> None:
        """ Check if it has any incoming node that is bridge """
        return any(node.is_bridge() for node in self.in_nodes)

    def get_variants_at(self, start:int, end:int=-1) -> List[seqvar.VariantRecord]:
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

    def get_cleavage_gain_from_downstream(self) -> List[seqvar.VariantRecord]:
        """ Get the variants that gains the cleavage by downstream nodes """
        cleavage_gain = []
        for node in self.out_nodes:
            if not node.variants:
                return []
            if node.variants[0].location.start != 0:
                return []
            if not cleavage_gain:
                cleavage_gain.append(node.variants[0].variant)
        return cleavage_gain

    def get_stop_lost_variants(self, stop_index:int) -> List[seqvar.VariantRecord]:
        """ Get stop lost variants """
        stop_lost = []
        stop_codon = FeatureLocation(start=stop_index, end=stop_index + 3)
        for variant in self.variants:
            if variant.variant.location.overlaps(stop_codon):
                stop_lost.append(variant.variant)
        return stop_lost

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

    def split_node(self, index:int, cleavage:bool=False, pre_cleave:bool=False,
            pop_collapse:bool=False) -> PVGNode:
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
        for variant in self.variants:
            if variant.location.start < index:
                if variant.location.end <= index:
                    left_variants.append(variant)
                else:
                    left_variants.append(variant[:index])
            if variant.location.end > index:
                right_variants.append(variant.shift(-index))
        self.seq = left_seq
        self.variants = left_variants

        new_node = PVGNode(
            seq=right_seq,
            reading_frame_index=self.reading_frame_index,
            variants=right_variants,
            orf=self.orf,
            was_bridge=self.was_bridge,
            pre_cleaved=pre_cleave,
            subgraph_id=self.subgraph_id,
            level=self.level
        )
        new_node.orf = self.orf

        if cleavage:
            new_node.cleavage = True

        if pop_collapse:
            self.cpop_collapsed = True
            new_node.npop_collapsed = True
            if new_node.is_bridge():
                self.was_bridge = True

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

        for variant in self.variants:
            if variant.location.start < i:
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))

        right_node = PVGNode(
            seq=right_seq,
            variants=right_variants,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            cpop_collapsed=self.cpop_collapsed
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

        for variant in self.variants:
            if variant.location.start < i:
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))

        left_node = PVGNode(
            seq=left_seq,
            variants=left_variants,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            npop_collapsed=self.npop_collapsed
        )

        self.seq = self.seq[i:]
        self.variants = right_variants

        return left_node

    def append_left(self, other:PVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq
        variants = copy.copy(other.variants)
        for variant in self.variants:
            variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants

    def append_right(self, other:PVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq
        for variant in other.variants:
            self.variants.append(variant.shift(len(self.seq.seq)))
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
            cleavage=self.cleavage,
            truncated=self.truncated,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            npop_collapsed=self.npop_collapsed,
            cpop_collapsed=self.cpop_collapsed
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

    def is_identical(self, other:PVGNode) -> bool:
        """ Checks if two nodes have the same sequence and same outbond node """
        return self.seq == other.seq and \
            self.out_nodes == other.out_nodes and \
            self.cleavage == other.cleavage and \
            self.reading_frame_index == other.reading_frame_index

    def is_less_mutated_than(self, other:PVGNode) -> bool:
        """ Checks if this node has less mutation than the other """
        if self.level != other.level:
            return self.level < other.level

        if len(self.variants) != len(other.variants):
            return len(self.variants) < len(other.variants)

        for x, y in zip(self.variants, other.variants):
            if x.location != y.location:
                return x.location < y.location

        for x, y in zip(self.variants, other.variants):
            # SNV is higher than RES
            if x.variant.type == 'SNV' and y.variant.type == 'RNAEditingSite':
                return True
            if x.variant.type == 'RNAEditingSite' and y.variant.type == 'SNV':
                return False
            # Otherwise don't really care about the order, just to make sure
            # the result is reproducible
            if x.variant.type != y.variant.type:
                return x.variant.type > y.variant.type

            if x.variant.alt != y.variant.type:
                return x.variant.alt > y.variant.type
        return True

    def transfer_in_nodes_to(self, other:PVGNode):
        """ transfer all in nodes of current node to the other """
        for in_node in copy.copy(self.in_nodes):
            if in_node in other.in_nodes:
                continue
            in_node.add_out_edge(other)

    def is_already_cleaved(self):
        """ Check if the node is already cleaved """
        return (self.cleavage and not self.pre_cleave) and \
            all(x.cleavage and not x.pre_cleave for x in self.out_nodes) and\
            all(all(y.cleavage and not y.pre_cleave for y in x.in_nodes)
                for x in self.out_nodes)
