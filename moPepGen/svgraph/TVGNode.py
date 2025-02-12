""" Module for TVGNode class """
from __future__ import annotations
from typing import List, Set, Tuple, Dict, Deque, TYPE_CHECKING
import copy
from collections import deque
import math
import uuid
from Bio.Seq import Seq
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen import seqvar, aa
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from .PVGNode import PVGNode
from .TVGEdge import TVGEdge

if TYPE_CHECKING:
    from moPepGen.svgraph import SubgraphTree

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
            level:int=0, was_bridge:bool=False, id:str=None):
        """ Constructor for TVGNode.

        Args:
            seq (DNASeqRecord): The sequence.
            variant (VariantRecord | None): The variant record or None for
                reference.
        """
        self.seq = seq
        self.in_edges:Set[TVGEdge] = set()
        self.out_edges:Set[TVGEdge] = set()
        self.variants:List[seqvar.VariantRecordWithCoordinate] = variants if variants else []
        self.branch = branch
        self.orf = orf or [None, None]
        self.reading_frame_index = reading_frame_index
        self.subgraph_id = subgraph_id
        self.global_variant = global_variant
        self.level = level
        self.was_bridge = was_bridge
        self.id = id or str(uuid.uuid4())

    def create_node(self, seq:DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            branch:bool=False, orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None, global_variant:seqvar.VariantRecord=None,
            level:int=None, was_bridge:bool=False):
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
            level=level or self.level,
            was_bridge=was_bridge
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
                new_start = max(variant.location.start, location.start)
                new_end = min(variant.location.end, location.end)
                new_start_offset = variant.location.start_offset \
                    if new_start == variant.location.start else 0
                new_end_offset = variant.location.end_offset \
                    if new_end == variant.location.end else 0
                new_loc = FeatureLocation(
                    start=max(variant.location.start, location.start),
                    end=min(variant.location.end, location.end),
                    seqname=variant.location.seqname,
                    reading_frame_index=variant.location.reading_frame_index,
                    start_offset=new_start_offset,
                    end_offset=new_end_offset
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

    def is_bridge(self) -> None:
        """ Check if this is a bridge node to another reading frame """
        for edge in self.out_edges:
            node = edge.out_node
            if node.reading_frame_index != self.reading_frame_index:
                return True
        if self.was_bridge:
            return True
        return False

    def get_edge_to(self, other:TVGNode) -> TVGEdge:
        """ Find the edge from this to the other node """
        for edge in self.out_edges:
            if edge.out_node is other:
                return edge
        raise ValueError('TVGEdge not found')

    def get_edge_from(self, other:TVGNode) -> TVGEdge:
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

    def _get_nth_rf_index(self, i:int) -> int:
        """ """
        if (i > 0 or i < -1) and not self.has_multiple_segments():
            raise ValueError('Node does not have multiple segments')

        locations = [loc.query for loc in self.seq.locations]
        locations += [v.location for v in self.variants
            if self.global_variant is None or v.variant != self.global_variant]

        locations.sort()

        return locations[i].reading_frame_index

    def get_first_rf_index(self) -> int:
        """ Get the first fragment's reading frame index """
        return self._get_nth_rf_index(0)

    def get_second_rf_index(self) -> int:
        """ Get the second fragment's reading frame index """
        return self._get_nth_rf_index(1)

    def get_last_rf_index(self) -> int:
        """ Get the last fragment's reading frame index """
        return self._get_nth_rf_index(-1)

    def _get_nth_subgraph_id(self, i) -> str:
        """ Get the nth fragment's subgraph ID """
        if (i > 0 or i < -1) and not self.has_multiple_segments():
            raise ValueError('Node does not have multiple segments')

        locations = [(loc.query, loc.ref.seqname) for loc in self.seq.locations]
        locations += [(v.location, v.location.seqname) for v in self.variants
            if self.global_variant is None or v.variant != self.global_variant]

        locations = sorted(locations, key=lambda x: x[0])

        return locations[i][1]

    def get_max_subgraph_id(self, subgraphs:SubgraphTree, last:bool=False) -> str:
        """ Get the max subgraph ID """
        max_subgraph_id = None
        max_level = -1
        if last:
            comp = lambda x,y: x >= y
        else:
            comp = lambda x,y: x > y
        for loc in self.seq.locations:
            subgraph_id = loc.ref.seqname
            level = subgraphs[subgraph_id].level
            if comp(level, max_level):
                max_level = level
                max_subgraph_id = subgraph_id
        for v in self.variants:
            subgraph_id = v.location.seqname
            level = subgraphs[subgraph_id].level
            if comp(level, max_level):
                max_subgraph_id = subgraph_id

        if not max_subgraph_id:
            max_subgraph_id = self.subgraph_id

        return max_subgraph_id

    def get_min_subgraph_id(self, subgraphs:SubgraphTree) -> str:
        """ Get the min subgraph ID """
        min_subgraph_id = None
        min_level = float('inf')
        for loc in self.seq.locations:
            subgraph_id = loc.ref.seqname
            level = subgraphs[subgraph_id].level
            if level < min_level:
                min_subgraph_id = subgraph_id
        for v in self.variants:
            subgraph_id = v.location.seqname
            level = subgraphs[subgraph_id].level
            if level < min_level:
                min_subgraph_id = subgraph_id

        if not min_subgraph_id:
            min_subgraph_id = self.subgraph_id

        return min_subgraph_id

    def get_first_subgraph_id(self) -> str:
        """ Get the first fragment's subgraph ID """
        return self._get_nth_subgraph_id(0)

    def get_last_subgraph_id(self) -> str:
        """ Get the last fragment's subgraph ID """
        return self._get_nth_subgraph_id(-1)

    def has_multiple_segments(self) -> bool:
        """ Whether the node has multiple segments, which is when the node
        is merged from several individual nodes. """
        n = len(self.seq.locations)
        n += len([v for v in self.variants if self.global_variant is not None
            and v.variant != self.global_variant])
        return n > 1

    def is_inbond_of(self, node:TVGEdge) -> bool:
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
        if len({x.variant for x in self.variants}) == 1 \
                and any(x.variant.is_real_fusion for x in self.variants):
            return True
        if self.global_variant is None:
            return not self.variants
        return not any(v.variant is not self.global_variant and
            v.variant != self.global_variant for v in self.variants)

    def is_stop_node(self) -> bool:
        """ check if it is a stop node """
        return self.seq.seq == '' and not self.out_edges

    def has_in_bridge(self) -> bool:
        """ check if it has any in node from different reading frame """
        return any(e.in_node.reading_frame_index != self.reading_frame_index
            for e in self.in_edges)

    def has_in_subgraph(self) -> bool:
        """ check if it has any in node from different reading frame """
        return any(e.in_node.subgraph_id != self.subgraph_id
            for e in self.in_edges)

    def is_bridge_to_subgraph(self) -> bool:
        """ check if it is a bridge node to a subgraph """
        return any(e.in_node.subgraph_id != self.subgraph_id for e in self.in_edges)

    def is_subgraph_bridge(self, out_node:TVGNode) -> bool:
        """ check if it is a subgraph bridge node """
        return self.subgraph_id != out_node.subgraph_id

    def is_subgraph_end(self) -> bool:
        """ check if is the end of a subgraph """
        return len(self.get_out_nodes()) > 0 \
            and all(x.level < self.level for x in self.get_out_nodes())

    def is_orf_bridge(self, out_node:TVGNode) -> bool:
        """ check if it is a orf bridge node """
        return self.reading_frame_index != out_node.reading_frame_index

    def is_out_orf_bridge(self) -> bool:
        """ check if it is an ORF out bridge node """
        return any(x.reading_frame_index != self.reading_frame_index
            for x in self.get_out_nodes())

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

    def is_exclusively_outbond_of(self, other:TVGNode) -> bool:
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
        variants = []
        for variant in self.variants:
            new_variant = seqvar.VariantRecordWithCoordinate(
                variant=variant.variant,
                location=copy.deepcopy(variant.location)
            )
            variants.append(new_variant)
        return TVGNode(
            seq=copy.deepcopy(self.seq),
            variants=variants,
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            global_variant=self.global_variant,
            level=self.level,
            was_bridge=self.was_bridge
        )

    def deepcopy(self, subgraph_id:str=None, level_increment:int=None) -> TVGNode:
        """ Create a deep copy of the node and all its downstream nodes.
        """
        new_node = self.copy()
        if level_increment is not None:
            new_node.level += level_increment

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
                        if new_out_node.seq:
                            for loc in new_out_node.seq.locations:
                                loc.ref.seqname = subgraph_id
                        for variant in new_out_node.variants:
                            variant.location.seqname = subgraph_id
                    if level_increment is not None:
                        new_out_node.level += level_increment
                    visited[source_out_node] = new_out_node
                new_edge = TVGEdge(target, new_out_node,
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
                end = min(i, variant.location.end)
                end_offset = variant.location.end_offset \
                    if end == variant.location.end else 0
                new_loc = FeatureLocation(
                    start=variant.location.start,
                    end=end,
                    seqname=variant.location.seqname,
                    reading_frame_index=variant.location.reading_frame_index,
                    start_offset=variant.location.start_offset,
                    end_offset=end_offset
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                left_variants.append(new_variant)
            if variant.location.end > i:
                start = max(variant.location.start, i)
                start_offset = variant.location.start_offset \
                    if start == variant.location.start else 0
                new_loc = FeatureLocation(
                    start=start,
                    end=variant.location.end,
                    seqname=variant.location.seqname,
                    reading_frame_index=variant.location.reading_frame_index,
                    start_offset=start_offset,
                    end_offset=variant.location.end_offset
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
                end = min(i, variant.location.end)
                end_offset = variant.location.end_offset \
                    if end == variant.location.end else 0
                new_loc = FeatureLocation(
                    start=variant.location.start,
                    end=end,
                    seqname=variant.location.seqname,
                    reading_frame_index=variant.location.reading_frame_index,
                    start_offset=variant.location.start_offset,
                    end_offset=end_offset
                )
                new_variant = seqvar.VariantRecordWithCoordinate(
                    variant=variant.variant, location=new_loc
                )
                left_variants.append(new_variant)
            if variant.location.end > i:
                start = max(variant.location.start, i)
                start_offset = variant.location.start_offset \
                    if start == variant.location.start else 0
                new_loc = FeatureLocation(
                    start=start,
                    end=variant.location.end,
                    seqname=variant.location.seqname,
                    reading_frame_index=variant.location.reading_frame_index,
                    start_offset=start_offset,
                    end_offset=variant.location.end_offset
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
            level=self.level,
            was_bridge=self.was_bridge
        )

        self.seq = self.seq[:i]
        self.variants = left_variants

        return right_node

    def append_left(self, other:TVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq

        variants = copy.copy(other.variants)
        for v_r in self.variants:
            v_l = None
            should_combine_variants = False

            for v_l in other.variants:
                if v_r.location.start != 0:
                    break
                if v_l.location.end != len(other.seq.seq):
                    continue
                if v_l.variant == v_r.variant:
                    should_combine_variants = False
                    break

            if should_combine_variants and v_l:
                v_l.location = FeatureLocation(
                    start=v_l.location.start,
                    end=v_r.location.end + len(other.seq.seq),
                    reading_frame_index=v_l.location.reading_frame_index,
                    start_offset=v_l.location.start_offset,
                    end_offset=v_r.location.end_offset
                )
            else:
                variants.append(v_r.shift(len(other.seq.seq)))
        self.variants = variants

    def append_right(self, other:TVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq

        for v_r in other.variants:
            v_l = None
            should_combine_variants = False
            for v_l in self.variants:
                if v_r.location.start != 0:
                    break
                if v_l.location.end != len(self.seq.seq):
                    continue
                if v_l.variant == v_r.variant \
                        and v_l.variant.is_real_fusion == v_r.variant.is_real_fusion:
                    should_combine_variants = True
                    break

            if should_combine_variants and v_l:
                v_l.location = FeatureLocation(
                    start=v_l.location.start,
                    end=v_r.location.end + len(self.seq.seq),
                    seqname=v_l.location.seqname,
                    reading_frame_index=v_l.location.reading_frame_index,
                    start_offset=v_l.location.start_offset,
                    end_offset=v_r.location.end_offset
                )
            else:
                self.variants.append(v_r.shift(len(self.seq.seq)))

        self.seq = new_seq

    def translate(self) -> PVGNode:
        """ translate to a PVGNode """
        if not self.out_edges:
            seq = self.seq[:len(self.seq) - len(self.seq) % 3].translate()
        else:
            seq = self.seq.translate()

        locations = []
        for loc in self.seq.locations:
            query_start = math.floor(loc.query.start/ 3)
            query_end = math.ceil(loc.query.end / 3)
            query_start_offset = loc.query.start - query_start * 3
            query_end_offset = query_end * 3 - loc.query.end
            query = FeatureLocation(
                start=query_start, end=query_end,
                reading_frame_index=loc.query.reading_frame_index,
                start_offset=query_start_offset,
                end_offset=query_end_offset
            )
            dna_query_codon_start = query_start * 3
            dna_ref_codon_start = loc.ref.start - (loc.query.start - dna_query_codon_start)
            ref_start = math.floor(dna_ref_codon_start / 3)
            dna_query_codon_end = query_end * 3
            dna_ref_codon_end = loc.ref.end + (dna_query_codon_end - loc.query.end)
            ref_end = ref_start + len(query)
            ref_start_offset = dna_ref_codon_start - ref_start * 3
            ref_end_offset = dna_ref_codon_end - ref_end * 3
            ref = FeatureLocation(
                start=ref_start, end=ref_end, seqname=loc.ref.seqname,
                start_offset=ref_start_offset, end_offset=ref_end_offset
            )
            locations.append(MatchedLocation(query=query, ref=ref))

        seq.__class__ = aa.AminoAcidSeqRecordWithCoordinates
        seq.locations = locations
        seq.orf = self.orf

        # translate the dna variant location to peptide coordinates.
        variants = [v.to_protein_coordinates() for v in self.variants]

        return PVGNode(
            seq=seq,
            variants=variants,
            orf=[None, None],
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id,
            level=self.level
        )

    def get_ith_variant_var_aa(self, i:int) -> Seq:
        """ Get the variant amino acid sequence of the ith variant """
        v = self.variants[i]
        loc = v.location
        lhs = loc.start - loc.start % 3
        rhs = loc.end + (3 - loc.end % 3) % 3
        rhs = min(rhs, len(self.seq.seq))
        seq = self.seq.seq[lhs:rhs]
        seq = seq[:len(seq) - len(seq) % 3]
        return seq.translate(to_stop=False)

    def get_ith_variant_ref_aa(self, i:int, tx_seq:Seq) -> Seq:
        """ Get the reference amino acid sequence of the ith variant. """
        v = self.variants[i]
        if v.variant.type == 'Insertion':
            left_offset = (v.location.start - 1) % 3
        else:
            left_offset = v.location.start % 3
        lhs = v.variant.location.start - left_offset

        seq = Seq('')
        if v.variant.type in ['Deletion', 'Substitution']:
            seq += tx_seq[v.variant.location.start:v.variant.location.end]
        else:
            seq += Seq(v.variant.ref)

        if left_offset > 0:
            seq = tx_seq[lhs:v.variant.location.start] + seq

        rhs = v.variant.location.end + (3 - len(seq) % 3) % 3
        rhs = min(rhs, len(tx_seq))

        if rhs > v.variant.location.end:
            seq += tx_seq[v.variant.location.end:rhs]
        seq = seq[:len(seq) - len(seq) % 3]
        return seq.translate(to_stop=False)

    def get_ref_sequence(self, tx_seq:Seq) -> Seq:
        """ Get the reference sequence """
        if not self.variants:
            return self.seq.seq

        seq = self.seq.seq[:self.variants[0].location.start]
        for i, variant in enumerate(self.variants):
            if variant.variant.type in ['Deletion', 'Substitution']:
                seq += tx_seq[variant.variant.location.start:variant.variant.location.end]
            elif variant.variant.type != 'Insertion':
                seq += variant.variant.ref
            if i + 1 >= len(self.variants):
                seq += self.seq.seq[variant.location.end:]
            else:
                next_var = self.variants[i + 1]
                seq += self.seq.seq[variant.location.end:next_var.location.start]

        return seq

    def check_stop_altering(self, tx_seq:Seq, cds_end:int=None):
        """ Checks if any variant is stop altering """
        if not self.variants:
            return

        if cds_end:
            stop_codon = FeatureLocation(start=cds_end, end=cds_end + 3)
            if not self.get_node_start_reference_index() <= cds_end < \
                    self.get_node_end_reference_index():
                return
            for variant in self.variants:
                if stop_codon.overlaps(variant.variant.location):
                    variant.is_stop_altering = True
            return

        if self.global_variant and self.global_variant.is_fusion():
            return

        for i,v in enumerate(self.variants):
            if v.variant.type in {'Fusion', 'circRNA'}:
                continue
            ref_aa = self.get_ith_variant_ref_aa(i, tx_seq)
            var_aa = self.get_ith_variant_var_aa(i)
            v.is_silent = v.variant.is_snv() and ref_aa == var_aa
            v.is_stop_altering = \
                (v.variant.is_snv() and ref_aa == '*' and var_aa != '*') \
                or (v.variant.is_insertion() and ref_aa =='*'
                    and not var_aa.startswith('*')) \
                or (v.variant.is_deletion() and '*' in ref_aa
                    and not (ref_aa.startswith('*') and var_aa.startswith('*'))) \
                or (v.variant.type in ['Substitution', 'MNV']
                    and '*' in ref_aa and ref_aa != var_aa)

    def is_less_mutated_than(self, other:TVGNode) -> bool:
        """ Checks if this node has less mutation than the other """
        self_vars = [v for v in self.variants if v.variant is not self.global_variant]
        other_vars = [v for v in other.variants if v.variant is not other.global_variant]
        if self.level != other.level:
            return self.level < other.level

        if len(self_vars) != len(other_vars):
            return len(self_vars) < len(other_vars)

        for x, y in zip(self_vars, other_vars):
            if x.location != y.location:
                return x.location < y.location

        for x, y in zip(self_vars, other_vars):
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

    def is_inframe_subgraph(self, start:TVGNode, end:TVGNode) -> bool:
        """ Check whether it is a inframe subgraph, i.e. its a subgraph node
        and its upstream and downstream are in the same reading frame and
        are also from the main graph. """
        if len(self.in_edges) != 1 or len(self.out_edges) != 1:
            return False

        upstream = self.get_in_nodes()[0]
        downstream = self.get_out_nodes()[0]

        return downstream.subgraph_id == start.subgraph_id \
            and upstream.subgraph_id == start.subgraph_id \
            and downstream.reading_frame_index == start.reading_frame_index \
            and upstream.reading_frame_index == start.reading_frame_index \
            and (upstream.seq.seq == ''
                or (upstream.seq.locations[0].ref > start.seq.locations[0].ref)) \
            and downstream.seq.locations[0].ref < end.seq.locations[0].ref

    def get_selenocysteine_positions(self, selenocysteines:List[FeatureLocation]
            ) -> List[int]:
        """ Find selenocysteine position from the sequence """
        positions = []
        for sec in selenocysteines:
            query_i = self.seq.get_query_index(sec.start)
            if query_i == -1:
                continue
            if query_i % 3 != 0:
                continue
            sec_query = FeatureLocation(start=query_i, end=query_i + 3)
            if any(v.location.overlaps(sec_query) for v in self.variants):
                continue
            positions.append(query_i)
        return positions
