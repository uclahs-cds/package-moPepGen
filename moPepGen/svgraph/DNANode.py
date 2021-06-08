""" Module for DNANode class """
from __future__ import annotations
from typing import List, Set, Tuple, Dict, Deque
from collections import deque
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen import seqvar, svgraph


class DNANode():
    """ Defines the nodes in the TranscriptVariantGraph

    Attributes:
        in_edges (Set[DNAEdge]): The inbonding edges.
        out_edges (Set[DNAEdge]): The outbonding edges.
        seq (DNASeqRecord): The sequence.
        variant (VariantRecord | None): The variant record or None for
            reference.
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            frameshifts:Set[seqvar.VariantRecord]=None):
        """ Constructor for DNANode.

        Args:
            seq (DNASeqRecord): The sequence.
            variant (VariantRecord | None): The variant record or None for
                reference.
        """
        self.seq = seq
        self.in_edges = set()
        self.out_edges = set()
        self.variants = variants if variants else []
        self.frameshifts = frameshifts if frameshifts else set()

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

    def get_edge_to(self, other:DNANode) -> svgraph.DNAEdge:
        """ Find the edge from this to the other node """
        for edge in self.out_edges:
            if edge.out_node is other:
                return edge
        raise ValueError('DNAEdge not found')

    def get_edge_from(self, other:DNANode) -> svgraph.DNAEdge:
        """ Find the edge from the other to this node. """
        for edge in self.in_edges:
            if edge.in_node is other:
                return edge
        raise ValueError('DNAEdge not found')

    def is_orphan(self) -> bool:
        """ Checks if the node is an orphan, that doesn't have any inbonding
        or outbonding edge. """
        return (not self.in_edges) and (not self.out_edges)

    def get_reference_next(self) -> DNANode:
        """ Get the next node of which the edge is reference (not variant
        or cleavage). """
        if not self.out_edges:
            return None
        for edge in self.out_edges:
            if edge.type in ['reference', 'variant_end']:
                return edge.out_node
        raise ValueError('No reference edge was found.')

    def get_reference_prev(self) -> DNANode:
        """ Get the previous node of which the edge is reference (not variant
        or cleavage) """
        for edge in self.in_edges:
            if edge.type == 'reference':
                return edge.in_node
        return None

    def deepcopy(self) -> DNANode:
        """ Create a deep copy of the node and all its downstream nodes. """

        new_node = self.__class__(
            seq=self.seq,
            variants=self.variants,
            frameshifts=self.frameshifts
        )

        queue:Deque[Tuple[DNANode, DNANode]] = deque([(self, new_node)])
        visited:Dict[DNANode, DNANode] = {}

        while queue:
            source, target = queue.pop()
            if source in visited:
                continue
            for edge in source.out_edges:
                source_out_node = edge.out_node
                if source_out_node in visited:
                    new_out_node = visited[source_out_node]
                else:
                    new_out_node = self.__class__(
                        seq=source_out_node.seq,
                        variants=source_out_node.variants,
                        frameshifts=source_out_node.frameshifts
                    )
                    visited[source_out_node] = new_out_node
                new_edge = svgraph.DNAEdge(target, new_out_node,
                    _type=edge.type)
                target.out_edges.add(new_edge)
                new_out_node.in_edges.add(new_edge)
                queue.appendleft((edge.out_node, new_out_node))

        return new_node

    def find_farthest_node_with_overlap(self) -> DNANode:
        r""" Find the farthest node, that within the range between the current
        node and it, there is at least one varint at any position of the
        reference sequence.

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
        if not self.get_reference_next().out_edges:
            return farthest
        queue = deque([edge.out_node for edge in self.out_edges])
        visited = set([self])
        while queue:
            cur:DNANode = queue.popleft()

            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                if not queue and cur is farthest:
                    # if the farthest has less than 5 neucleotides, continue
                    # searching, because it's likely not able to keep one amino
                    # acid after translation.
                    if len(cur.seq) < 5:
                        for edge in cur.out_edges:
                            queue.append(edge.out_node)
                continue

            if farthest is None and not cur.variants:
                farthest = cur
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
