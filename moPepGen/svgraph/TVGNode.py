""" Module for TVGNode class """
from __future__ import annotations
from typing import List, Set, Tuple, Dict, Deque, TYPE_CHECKING
import copy
from collections import deque
from functools import reduce
import math
import re
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen import seqvar, svgraph, aa
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation


if TYPE_CHECKING:
    from Bio.Seq import Seq

class TVGNode():
    """ Defines the nodes in the TranscriptVariantGraph

    Attributes:
        in_edges (Set[TVGEdge]): The inbonding edges.
        out_edges (Set[TVGEdge]): The outbonding edges.
        seq (DNASeqRecord): The sequence.
        variant (VariantRecord | None): The variant record or None for
            reference.
        branch (bool)
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            branch:bool=False, orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None, global_variant:seqvar.VariantRecord=None,
            level:int=0):
        """ Constructor for TVGNode.

        Args:
            seq (DNASeqRecord): The sequence.
            variant (VariantRecord | None): The variant record or None for
                reference.
        """
        self.seq = seq
        self.in_edges:Set[svgraph.TVGEdge] = set()
        self.out_edges:Set[svgraph.TVGEdge] = set()
        self.variants:List[seqvar.VariantRecordWithCoordinate] = variants if variants else []
        self.branch = branch
        self.orf = orf or [None, None]
        self.reading_frame_index = reading_frame_index
        self.subgraph_id = subgraph_id
        self.global_variant = global_variant
        self.level = level

    def create_node(self, seq:DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            branch:bool=False, orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None, global_variant:seqvar.VariantRecord=None,
            level:int=None):
        """ Constructor for TVGNode.

        Args:
            seq (DNASeqRecord): The sequence.
            variant (VariantRecord | None): The variant record or None for
                reference.
        """
        return self.__class__(
            seq=seq,
            variants=variants,
            branch=branch,
            orf=orf,
            reading_frame_index=reading_frame_index,
            subgraph_id=subgraph_id,
            global_variant=global_variant or self.global_variant,
            level=level or self.level
        )

    def __hash__(self):
        """ hash """
        if self.seq:
            seq = self.seq.seq
            locations = [(it.ref.start, it.ref.end) for it in \
                self.seq.locations]
        else:
            seq = None
            locations = []
        return hash((str(seq), *locations))

    def __getitem__(self, index) -> TVGNode:
        """ get item """
        start, stop, _ = index.indices(len(self.seq))
        location = FeatureLocation(start=start, end=stop)
        seq = self.seq.__getitem__(index)
        variants = []

        for variant in self.variants:
            if variant.location.overlaps(location):
                new_loc = FeatureLocation(
                    start=max(variant.location.start, location.start),
                    end=min(variant.location.end, location.end)
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                new_variant = new_variant.shift(-start)
                variants.append(new_variant)
        return self.create_node(
            seq=seq,
            variants=variants,
            orf=self.orf,
            branch=self.branch,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            level=self.level
        )

    def get_edge_to(self, other:TVGNode) -> svgraph.TVGEdge:
        """ Find the edge from this to the other node """
        for edge in self.out_edges:
            if edge.out_node is other:
                return edge
        raise ValueError('TVGEdge not found')

    def get_edge_from(self, other:TVGNode) -> svgraph.TVGEdge:
        """ Find the edge from the other to this node. """
        for edge in self.in_edges:
            if edge.in_node is other:
                return edge
        raise ValueError('TVGEdge not found')

    def get_out_nodes(self) -> List[TVGNode]:
        """ get outbonding nodes as a list """
        return [x.out_node for x in self.out_edges]

    def get_in_nodes(self) -> List[TVGNode]:
        """ get incoming nodes as a list """
        return [x.in_node for x in self.in_edges]

    def is_inbond_of(self, node:svgraph.TVGEdge) -> bool:
        """ Checks if this node is a inbound node of the other node. """
        for edge in self.out_edges:
            if node is edge.out_node:
                return True
        return False

    def is_orphan(self) -> bool:
        """ Checks if the node is an orphan, that doesn't have any inbonding
        or outbonding edge. """
        return (not self.in_edges) and (not self.out_edges)

    def is_reference(self) -> bool:
        """ check if it is reference (no variants) """
        if len(self.variants) == 1 and self.variants[0].variant.is_real_fusion:
            return True
        if self.global_variant is None:
            return not self.variants
        return not any(v.variant is not self.global_variant for v in self.variants)

    def is_stop_node(self) -> bool:
        """ check if it is a stop node """
        return self.seq.seq == '' and not self.out_edges

    def has_in_bridge(self) -> bool:
        """ check if it has any in node from different reading frame """
        return any(e.in_node.reading_frame_index != self.reading_frame_index
            for e in self.in_edges)

    def is_bridge_to_subgraph(self) -> bool:
        """ check if it is a bridge node to a subgraph """
        return any(e.in_node.subgraph_id != self.subgraph_id for e in self.in_edges)

    def is_subgraph_bridge(self, out_node:TVGNode) -> bool:
        """ check if it is a subgraph bridge node """
        return self.subgraph_id != out_node.subgraph_id

    def is_orf_bridge(self, out_node:TVGNode) -> bool:
        """ check if it is a orf bridge node """
        return self.reading_frame_index != out_node.reading_frame_index

    def has_bridge_from_reading_frame(self, other_reading_frame_index:int):
        """ Check if it has a in bridge node from another reading frame """
        return any(x.reading_frame_index == other_reading_frame_index\
            for x in self.get_in_bridges())

    def has_ref_position(self, i:int) -> bool:
        """ Check if the node has a given reference position """
        return any(i in x.ref for x in self.seq.locations)

    def get_in_bridges(self) -> List[TVGNode]:
        """ Get all in bridge nodes """
        in_bridges = []
        for edge in self.in_edges:
            in_node = edge.in_node
            if any(v.variant.is_frameshifting() for v in in_node.variants):
                in_bridges.append(in_node)
        return in_bridges

    def is_exclusively_outbond_of(self, other:svgraph.TVGNode) -> bool:
        """ Checks if the node is exclusively outbond with the other.
        Exclusive binding means the upstream node has only 1 outbond edge,
        and the downstream node has only inbond edge."""
        if not self.is_inbond_of(other):
            return False
        return len(self.out_edges) == 1 and len(other.in_edges) == 1

    def has_exclusive_outbond_node(self) -> bool:
        """ The given node has exclusive outbond node """
        return len(self.out_edges) == 1 and \
            len(self.get_out_nodes()[0].in_edges) == 1

    def has_exclusive_inbond_node(self) -> bool:
        """ The given node has exclusive inbond node """
        return len(self.in_edges) == 1 and \
            len(self.get_in_nodes()[0].out_edges) == 1

    def get_reference_next(self) -> TVGNode:
        """ Get the next node of which the edge is reference (not variant
        or cleavage). """
        if not self.out_edges:
            return None
        if len(self.out_edges) == 1:
            return list(self.out_edges)[0].out_node
        ref_node = None
        for edge in self.out_edges:
            if edge.type in ['reference', 'variant_end']:
                if ref_node is None:
                    ref_node = edge.out_node
                elif not ref_node.branch:
                    ref_node = edge.out_node
        if not ref_node:
            raise ValueError('No reference edge was found.')
        return ref_node

    def get_reference_prev(self) -> TVGNode:
        """ Get the previous node of which the edge is reference (not variant
        or cleavage) """
        for edge in self.in_edges:
            if edge.type == 'reference':
                return edge.in_node
        return None

    def get_node_start_reference_index(self) -> int:
        """ Get the referece index of node start """
        if self.seq.locations and self.variants:
            return min(
                self.seq.locations[0].ref.start,
                self.variants[0].variant.location.start
            )

        if self.seq.locations:
            return self.variants[0].ref.start

        return self.variants[0].variant.location.start

    def get_node_end_reference_index(self) -> int:
        """ Get the reference index of node end """
        if self.seq.locations and self.variants:
            return max(
                self.seq.locations[-1].ref.end,
                self.variants[-1].variant.location.end
            )

        if self.seq.locations:
            return self.variants[-1].ref.end

        return self.variants[-1].variant.location.end

    def copy(self) -> TVGNode:
        """ Create a copy of the node """
        return TVGNode(
            seq=copy.copy(self.seq),
            variants=copy.copy(self.variants),
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            global_variant=self.global_variant,
            level=self.level
        )

    def deepcopy(self, subgraph_id:str=None, level_increment:int=None) -> TVGNode:
        """ Create a deep copy of the node and all its downstream nodes.
        """
        new_node = self.copy()

        queue:Deque[Tuple[TVGNode, TVGNode]] = deque([(self, new_node)])
        visited:Dict[TVGNode, TVGNode] = {}

        while queue:
            source, target = queue.pop()
            for edge in source.out_edges:
                source_out_node = edge.out_node
                visited_this = False
                if source_out_node in visited:
                    new_out_node = visited[source_out_node]
                    visited_this = True
                else:
                    new_out_node = source_out_node.copy()
                    if subgraph_id:
                        new_out_node.subgraph_id = subgraph_id
                    if level_increment is not None:
                        new_out_node.level += level_increment
                    visited[source_out_node] = new_out_node
                new_edge = svgraph.TVGEdge(target, new_out_node,
                    _type=edge.type)
                target.out_edges.add(new_edge)
                new_out_node.in_edges.add(new_edge)
                if not visited_this and edge.out_node.out_edges:
                    queue.appendleft((edge.out_node, new_out_node))

        return new_node

    def stringify(self, k:int=None) -> None:
        """ Get a str representation of a subgraph """
        if not k:
            k = float("inf")

        _type = 'alt' if self.variants else 'ref'
        if self.seq:
            n = len(self.seq.seq)
            seq_str = f'{str(self.seq.seq[:k])}...' if n >= k else self.seq.seq
            node = f'{seq_str}|{_type}'
        else:
            node = 'root'

        if not self.out_edges:
            return {node: '$'}

        downstream = {}
        for edge in self.out_edges:
            downstream.update(edge.out_node.stringify(k))

        return {node: downstream}

    def get_orf_start(self, i:int=0) -> int:
        """ Get the ORF start position given the start codon is found at
        position i of the sequence. """
        if self.seq.locations:
            for loc in self.seq.locations:
                if i < loc.query.start:
                    return int((3 - (loc.query.start - i) % 3) % 3 + loc.ref.start)
                if loc.query.start <= i < loc.query.end:
                    return int(i - loc.query.start + loc.ref.start)
        out_node = self.get_reference_next()
        if out_node.seq.locations:
            loc = out_node.seq.locations[0]
            return int(((len(self.seq.seq) + loc.query.start - i) % 3) % 3 + loc.ref.start)
        raise ValueError('Can not find ORF')

    def truncate_left(self, i:int) -> TVGNode:
        """ Truncate the left i nucleotides off. A new node with the left part
        of the sequences and variants associated is returned. The self node
        is updated with only the right part of the sequence and variants. """
        left_seq = self.seq[:i]
        left_variants = []
        right_variants = []

        for variant in self.variants:
            if variant.location.start < i:
                new_loc = FeatureLocation(
                    start=variant.location.start,
                    end=min(i, variant.location.end)
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                left_variants.append(new_variant)
            if variant.location.end > i:
                new_loc = FeatureLocation(
                    start=max(variant.location.start, i),
                    end=variant.location.end
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                new_variant = new_variant.shift(-i)
                right_variants.append(new_variant)

        left_node = self.create_node(
            seq=left_seq,
            variants=left_variants,
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            global_variant=self.global_variant,
            level=self.level
        )

        self.seq = self.seq[i:]
        self.variants = right_variants
        self.branch = False

        return left_node

    def truncate_right(self, i:int) -> TVGNode:
        """ Truncate the right i nucleotides off. A new node with the right
        part of the sequences and variants associated is returned. The self
        node is updated with only the left part of the sequence and variants.
        """
        right_seq = self.seq[i:]
        left_variants = []
        right_variants = []

        for variant in self.variants:
            if variant.location.start < i:
                new_loc = FeatureLocation(
                    start=variant.location.start,
                    end=min(i, variant.location.end)
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                left_variants.append(new_variant)
            if variant.location.end > i:
                new_loc = FeatureLocation(
                    start=max(variant.location.start, i),
                    end=variant.location.end
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                new_variant = new_variant.shift(-i)
                right_variants.append(new_variant)

        right_node = self.create_node(
            seq=right_seq,
            variants=right_variants,
            branch=False,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            level=self.level
        )

        self.seq = self.seq[:i]
        self.variants = left_variants

        return right_node

    def append_left(self, other:TVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq

        variants = copy.copy(other.variants)
        for variant in self.variants:
            should_combine_variants = variants and \
                variants[-1].location.end == variant.location.start and \
                variants[-1] == variant.variant
            if should_combine_variants:
                variant.location = FeatureLocation(
                    start=other.variants[-1].location.start,
                    end=variant.location.end
                )
            else:
                variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants

    def append_right(self, other:TVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq

        for variant in other.variants:
            should_combine_variants = self.variants and \
                self.variants[-1].location.end == variant.location.start and \
                self.variants[-1].variant == variant.variant
            if should_combine_variants:
                self.variants[-1].location = FeatureLocation(
                    start=self.variants[-1].location.start,
                    end=variant.location.end
                )
            else:
                self.variants.append(variant.shift(len(self.seq.seq)))

        self.seq = new_seq

    def translate(self) -> svgraph.PVGNode:
        """ translate to a PVGNode """
        seq = self.seq.translate()

        locations = []
        for loc in self.seq.locations:
            if len(loc.query) < 3:
                continue
            query_start = math.ceil(loc.query.start / 3)
            query_end = math.floor(loc.query.end / 3)
            query = FeatureLocation(start=query_start, end=query_end)
            dna_query_codon_start = query_start * 3
            dna_ref_codon_start = loc.ref.start + dna_query_codon_start - query.start
            ref_start = math.ceil((dna_ref_codon_start - self.reading_frame_index) / 3)
            ref_end = ref_start + len(query)
            ref = FeatureLocation(start=ref_start, end=ref_end)
            locations.append(MatchedLocation(query=query, ref=ref))

        seq.__class__ = aa.AminoAcidSeqRecordWithCoordinates
        seq.locations = locations
        seq.orf = self.orf

        # translate the dna variant location to peptide coordinates.
        variants = [v.to_protein_coordinates() for v in self.variants]

        return svgraph.PVGNode(
            seq=seq,
            variants=variants,
            orf=[None, None],
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            level=self.level
        )

    def get_ref_sequence(self) -> Seq:
        """ Get the reference sequence """
        if not self.variants:
            return self.seq.seq

        seq = self.seq.seq[:self.variants[0].location.start]
        for i, variant in enumerate(self.variants):
            seq += variant.variant.ref
            if i + 1 >= len(self.variants):
                seq += self.seq.seq[variant.location.end:]
            else:
                next_var = self.variants[i + 1]
                seq += self.seq.seq[variant.location.end:next_var.location.start]

        return seq

    def check_stop_altering(self, cds_end:int=None):
        """ Checks whether each variant is stop lost """
        if not self.variants:
            return

        if cds_end:
            stop_codon = FeatureLocation(start=cds_end, end=cds_end + 3)
            if not self.get_node_start_reference_index() < cds_end < \
                    self.get_node_end_reference_index():
                return
            for variant in self.variants:
                if stop_codon.overlaps(variant.variant.location):
                    variant.is_stop_altering = True
            return

        ref_seq = self.get_ref_sequence()
        ref_seq = ref_seq[:math.floor(len(ref_seq)/3)*3]
        ref_aa = ref_seq.translate(to_stop=False)
        stop_positions = [x.start() for x in re.finditer(r'\*', str(ref_aa))]
        if not stop_positions:
            return

        shifted = reduce(
            lambda x,y: x+y,
            [len(x.variant.alt) - len(x.variant.ref) for x in self.variants]
        )

        for stop in stop_positions:
            pos = stop * 3 + shifted
            stop_codon = FeatureLocation(start=pos, end=pos + 3)
            for variant in self.variants:
                if stop_codon.overlaps(variant.location):
                    variant.is_stop_altering = True
