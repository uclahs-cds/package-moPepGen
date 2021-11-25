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
        frameshifts (seqvar.VariatnRecord)
        branch (bool)
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            frameshifts:Set[seqvar.VariantRecord]=None, branch:bool=False,
            orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None):
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
        self.frameshifts = frameshifts or set()
        self.branch = branch
        self.orf = orf or [None, None]
        self.reading_frame_index = reading_frame_index
        self.subgraph_id = subgraph_id

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
        frameshifts = copy.copy(self.frameshifts)

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
            elif variant.location.start >= stop and variant.varint in frameshifts:
                frameshifts.remove(variant.variant)
        return TVGNode(
            seq=seq,
            variants=variants,
            frameshifts=frameshifts,
            orf=self.orf,
            branch=self.branch,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id
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
        return not self.variants

    def has_in_bridge(self) -> bool:
        """ check if it has any in node from different reading frame """
        return any(e.in_node.reading_frame_index != self.reading_frame_index
            for e in self.in_edges)

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
            seq=self.seq,
            variants=copy.copy(self.variants),
            frameshifts=copy.copy(self.frameshifts),
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id
        )

    def deepcopy(self) -> TVGNode:
        """ Create a deep copy of the node and all its downstream nodes.

        Args:
            propagate_frameshifts (bool):
        """
        new_node = self.__class__(
            seq=self.seq,
            variants=copy.copy(self.variants),
            frameshifts=copy.copy(self.frameshifts),
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id
        )

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
                    new_out_node = source_out_node.__class__(
                        seq=source_out_node.seq,
                        variants=copy.copy(source_out_node.variants),
                        frameshifts=copy.copy(source_out_node.frameshifts),
                        branch=source_out_node.branch,
                        orf=source_out_node.orf,
                        reading_frame_index=source_out_node.reading_frame_index,
                        subgraph_id=source_out_node.subgraph_id
                    )
                    visited[source_out_node] = new_out_node
                new_edge = svgraph.TVGEdge(target, new_out_node,
                    _type=edge.type)
                target.out_edges.add(new_edge)
                new_out_node.in_edges.add(new_edge)
                if not visited_this and edge.out_node.out_edges:
                    queue.appendleft((edge.out_node, new_out_node))

        return new_node

    def find_farthest_node_with_overlap(self, min_size:int=6) -> TVGNode:
        r""" Find the farthest node, that within the range between the current
        node and it, there is at least one varint at any position of the
        reference sequence. If the farthest node found has an exclusive single
        out node, it extends to. For circular graph, this extension won't
        continue if the exclusive single node is root.

        For example, in a graph like below, the node ATGG's farthest node
        with overlap would be node 'CCCT'
                 T--
                /   \
            ATGG-TCT-G-CCCT
                    \ /
                     A
        """
         # find the range of overlaps
        farthest = None
        if not self.get_reference_next():
            return None
        if not self.get_reference_next().out_edges:
            return self.get_reference_next()
        queue = deque([edge.out_node for edge in self.out_edges])
        visited = {self}
        while queue:
            cur:TVGNode = queue.popleft()
            if cur is None:
                continue

            if cur.reading_frame_index != self.reading_frame_index:
                continue

            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                if cur is farthest and cur is not self:
                    if cur.out_edges and cur.get_reference_next().has_in_bridge():
                        for edge in cur.out_edges:
                            queue.append(edge.out_node)
                    # if the farthest has less than 6 neucleotides, continue
                    # searching, because it's likely not able to keep one amino
                    # acid after translation.
                    if len(cur.seq) < min_size:
                        for edge in cur.out_edges:
                            queue.append(edge.out_node)
                continue

            if farthest is None and not cur.variants:
                farthest = cur
                queue.append(cur)
                continue

            if cur.variants:
                next_ref = cur.get_reference_next()
                queue.append(next_ref)
                continue

            if cur.seq.locations[0] < farthest.seq.locations[0]:
                for edge in cur.out_edges:
                    queue.append(edge.out_node)
                continue

            if cur.seq.locations[0] > farthest.seq.locations[0]:
                farthest, cur = cur, farthest
                queue.append(cur)
                visited.remove(cur)
                continue
        return farthest

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
        left_frameshifts = copy.copy(self.frameshifts)

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
                if variant.location.start >= i and \
                        variant.variant.is_frameshifting():
                    left_frameshifts.remove(variant.variant)

        left_node = TVGNode(
            seq=left_seq,
            variants=left_variants,
            frameshifts=left_frameshifts,
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id
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
        right_frameshifts = copy.copy(self.frameshifts)

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
                if variant.location.start >= i and \
                        variant.variant.is_frameshifting():
                    self.frameshifts.discard(variant.variant)

        right_node = TVGNode(
            seq=right_seq,
            variants=right_variants,
            frameshifts=right_frameshifts,
            branch=False,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id
        )

        self.seq = self.seq[:i]
        self.variants = left_variants

        return right_node

    def append_left(self, other:TVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq

        frameshifts = copy.copy(other.frameshifts)
        frameshifts.update(self.variants)
        self.frameshifts = frameshifts

        variants = copy.copy(other.variants)
        for variant in self.variants:
            variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants

    def append_right(self, other:TVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq

        for variant in other.variants:
            self.variants.append(variant.shift(len(self.seq.seq)))
            if variant.variant.is_frameshifting():
                self.frameshifts.add(variant)

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
            frameshifts=copy.copy(self.frameshifts),
            orf=[None, None],
            reading_frame_index=self.reading_frame_index
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
