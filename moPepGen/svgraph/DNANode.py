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
        frameshifts (seqvar.VariatnRecord)
        branch (bool)
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            frameshifts:Set[seqvar.VariantRecord]=None, branch:bool=False):
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
        self.branch = branch

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

    def is_inbond_of(self, node:svgraph.DNAEdge) -> bool:
        """ Checks if this node is a inbound node of the other node. """
        for edge in self.out_edges:
            if node is edge.out_node:
                return True
        return False

    def is_orphan(self) -> bool:
        """ Checks if the node is an orphan, that doesn't have any inbonding
        or outbonding edge. """
        return (not self.in_edges) and (not self.out_edges)

    def needs_branch_out(self, branch_out_size:int=100) -> bool:
        """ Checks if the node needs to branch out """
        variant = self.variants[0].variant
        location = self.variants[0].location
        variant_len = location.end - location.start
        return variant.is_frameshifting() or variant_len >= branch_out_size

    def next_node_to_branch_out(self, to_node:svgraph.DNANode,
            branch_out_size:int=100) -> svgraph.DNANode:
        """ Returns the next node that needs to be branched out between two
        reference nodes. Returns None if not found.

        Args:
            to_node (svgraph.DNANode): The to node
            branch_out_size (int): The size limit of the variant to forch
                creating a branch in th graph.
        """
        queue:Deque[svgraph.DNANode] = deque()
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

    def deepcopy(self, propagate_frameshifts:bool=True) -> DNANode:
        """ Create a deep copy of the node and all its downstream nodes.

        Args:
            propagate_frameshifts (bool):
        """
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
                source_out_node:svgraph.DNANode = edge.out_node
                if source_out_node in visited:
                    new_out_node = visited[source_out_node]
                else:
                    frameshifts = source_out_node.frameshifts
                    if propagate_frameshifts:
                        frameshifts.update(self.frameshifts)
                    new_out_node = self.__class__(
                        seq=source_out_node.seq,
                        variants=source_out_node.variants,
                        frameshifts=frameshifts
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

            # if this is the start of a new branch, don't count it.
            if cur.branch:
                continue

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
