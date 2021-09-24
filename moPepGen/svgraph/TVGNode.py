""" Module for TVGNode class """
from __future__ import annotations
from typing import List, Set, Tuple, Dict, Deque
import copy
from collections import deque
import math
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen import seqvar, svgraph, aa
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation


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
            orf:List[int]=None):
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
        self.frameshifts = frameshifts if frameshifts else set()
        self.branch = branch
        self.orf = orf or [None, None]

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
                variants.append(variant.shift(-start))
            elif variant.location.start >= stop and variant.varint in frameshifts:
                frameshifts.remove(variant.variant)
        return TVGNode(
            seq=seq,
            variants=variants,
            frameshifts=frameshifts,
            orf=self.orf,
            branch=self.branch
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

    def should_skip_frameshift(self, dist:int=60, n:int=3) -> bool:
        """ checks if the node is frameshifting and should be skipped """
        my_frameshifts:List[seqvar.VariantRecord] = []
        for v in self.variants:
            if v.variant.is_frameshifting():
                my_frameshifts.append(v.variant)
        if not my_frameshifts:
            return False

        if len(self.frameshifts) > n:
            return True

        my_frameshifts.sort()
        this_frameshift = my_frameshifts[0]

        upstream_frameshifts:List[seqvar.VariantRecord] = []
        for v in self.frameshifts:
            if v not in my_frameshifts:
                upstream_frameshifts.append(v)
        if not upstream_frameshifts:
            return False
        last_frameshift = upstream_frameshifts[-1]

        n_shifted = sum([it.frames_shifted() for it in upstream_frameshifts])
        if n_shifted % 3 != 0:
            return False

        diff = this_frameshift.location.start - last_frameshift.location.end
        return diff >= dist

    def is_exclusively_outbond_of(self, other:svgraph.TVGNode) -> bool:
        """ Checks if the node is exclusively outbond with the other.
        Exclusive binding means the upstream node has only 1 outbond edge,
        and the downstream node has only inbond edge."""
        if not self.is_inbond_of(other):
            return False
        return len(self.out_edges) == 1 and len(other.in_edges) == 1

    def needs_branch_out(self, branch_out_size:int=100) -> bool:
        """ Checks if the node needs to branch out """
        for variant in self.variants:
            if variant.variant.is_frameshifting():
                return True
        for variant in self.variants:
            location = variant.location
            variant_len = location.end - location.start
            if variant_len >= branch_out_size:
                return True
        # For variants located in the end of the sequence.
        if not self.out_edges:
            return True
        return False

    def next_node_to_branch_out(self, to_node:svgraph.TVGNode,
            branch_out_size:int=100) -> svgraph.TVGNode:
        """ Returns the next node that needs to be branched out between two
        reference nodes. Returns None if not found.

        Args:
            to_node (svgraph.TVGNode): The to node
            branch_out_size (int): The size limit of the variant to forch
                creating a branch in th graph.
        """
        queue:Deque[svgraph.TVGNode] = deque()
        for edge in self.out_edges:
            queue.appendleft(edge.out_node)

        while queue:
            cur = queue.pop()
            if cur is to_node:
                continue
            if cur.branch:
                continue
            if not cur.variants:
                for edge in cur.out_edges:
                    queue.appendleft(edge.out_node)
                continue
            if cur.needs_branch_out(branch_out_size):
                return cur
            for edge in cur.out_edges:
                queue.appendleft(edge.out_node)

        return None

    def get_reference_next(self) -> TVGNode:
        """ Get the next node of which the edge is reference (not variant
        or cleavage). """
        if not self.out_edges:
            return None
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

    def copy(self) -> TVGNode:
        """ Create a copy of the node """
        return TVGNode(
            seq=self.seq,
            variants=copy.copy(self.variants),
            frameshifts=copy.copy(self.frameshifts),
            branch=self.branch,
            orf=self.orf
        )

    def deepcopy(self, propagate_frameshifts:bool=True) -> TVGNode:
        """ Create a deep copy of the node and all its downstream nodes.

        Args:
            propagate_frameshifts (bool):
        """
        new_node = self.__class__(
            seq=self.seq,
            variants=copy.copy(self.variants),
            frameshifts=copy.copy(self.frameshifts)
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
                    frameshifts = copy.copy(source_out_node.frameshifts)
                    if propagate_frameshifts:
                        frameshifts.update(target.frameshifts)
                    new_out_node = self.__class__(
                        seq=source_out_node.seq,
                        variants=copy.copy(source_out_node.variants),
                        frameshifts=frameshifts
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

            # if this is the start of a new branch, don't count it.
            if cur.branch:
                continue

            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                if cur is farthest and cur is not self:
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
            for loc in out_node.seq.locations:
                if i < loc.query.start:
                    return int((3 - (loc.query.start - i) % 3) % 3 + loc.ref.start)
                if loc.query.start <= i < loc.query.end:
                    return int(i - loc.query.start + loc.ref.start)
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
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))
                if variant.location.start >= i and \
                        variant.variant.is_frameshifting():
                    left_frameshifts.remove(variant.variant)

        left_node = TVGNode(
            seq=left_seq,
            variants=left_variants,
            frameshifts=left_frameshifts,
            branch=self.branch
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
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))
                if variant.location.start >= i and \
                        variant.variant.is_frameshifting():
                    self.frameshifts.remove(variant.variant)

        right_node = TVGNode(
            seq=right_seq,
            variants=right_variants,
            frameshifts=right_frameshifts,
            branch=False
        )

        self.seq = self.seq[:i]
        self.variants = left_variants

        return right_node

    def append_left(self, other:TVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq

        variants = copy.copy(other.variants)
        for variant in self.variants:
            variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants

        self.frameshifts.update(other.frameshifts)

    def append_right(self, other:TVGNode) -> None:
        """ Combine the other node the the right. """
        self.seq = self.seq + other.seq

        for variant in other.variants:
            self.variants.append(variant.shift(len(self.seq.seq)))

        self.frameshifts.update(other.frameshifts)

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
            ref_start = math.ceil((dna_ref_codon_start - self.orf[0]) / 3)
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
            orf=[None, None]
        )
