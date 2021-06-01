""" Module for DNANode class """
from __future__ import annotations
from typing import List, Set, Tuple
from collections import deque
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord
from moPepGen import svgraph


class DNANode():
    """ Defines the nodes in the TranscriptVariantGraph
    
    Attributes:
        in_edges (Set[DNAEdge]): The inbonding edges.
        out_edges (Set[DNAEdge]): The outbonding edges.
        seq (DNASeqRecord): The sequence.
        variant (VEPVariantRecord | None): The variant record or None for
            reference.
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates,
            variants:List[svgraph.VariantRecordWithCoordinate]=None,
            frameshifts:List[VEPVariantRecord]=None):
        """ Constructor for DNANode.
        
        Args:
            seq (DNASeqRecord): The sequence.
            variant (VEPVariantRecord | None): The variant record or None for
                reference.
        """
        self.seq = seq
        self.in_edges = set()
        self.out_edges = set()
        self.variants = variants if variants else []
        self.framshifts = frameshifts if frameshifts else []
    
    def __hash__(self):
        """ hash """
        locations = [(it.ref.start, it.ref.end) for it in self.seq.locations]
        return hash((str(self.seq.seq), *locations))
    
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
        for edge in self.out_edges:
            if edge.type in ['reference', 'variant_end']:
                return edge.out_node
        return None
    
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
            frameshifts=self.framshifts
        )

        queue:deque[Tuple[DNANode, DNANode]] = deque([(self, new_node)])

        while queue:
            source, target = queue.pop()
            for edge in source.out_edges:
                new_out_node = self.__class__(
                    seq=edge.out_node.seq,
                    variants=edge.out_node.variants,
                    frameshifts=edge.out_node.framshifts
                )
                new_edge = svgraph.DNAEdge(target, new_out_node, type=edge.type)
                target.out_edges.add(new_edge)
                new_out_node.in_edges.add(new_edge)
                queue.appendleft((edge.out_node, new_out_node))

        return new_node
    
    def find_farthest_node_with_overlap(self) -> DNANode:
        """ Find the farthest node, that within the range between the current
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
        queue = deque([edge.out_node for edge in self.out_edges])
        visited = set([self])
        while queue:
            cur:DNANode = queue.popleft()
            
            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
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
            
            if cur is farthest:
                continue
        return farthest