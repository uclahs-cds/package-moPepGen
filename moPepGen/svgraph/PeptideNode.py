""""""
from __future__ import annotations
import copy
from typing import List, Set
from moPepGen import vep
from moPepGen import aa
from moPepGen.svgraph.VariantRecordWithCoordinate import VariantRecordWithCoordinate


class PeptideNode():
    """ The PeptideNode class is used in the svgraph.PeptideVariantGraph, that
    stores the sequence in a Node.
    
    Attributes:
        seq (aa.AminoAcidSeqRecord): The amino acid sequence.
        variants (List[VariantRecordWithCoordinate]): The variant records
            carried in the sequence.
        frameshifts (Set[VEPVariantRecord]): Frameshifting variants.
        in_nodes (Set[PeptideNode]): Inbound nodes.
        ou_nodes (Set[PeptideNode]): Outbound nodes
    """
    def __init__(self, seq:aa.AminoAcidSeqRecord,
            variants:List[VariantRecordWithCoordinate]=None,
            in_nodes:Set[PeptideNode]=None,
            out_nodes:Set[PeptideNode]=None,
            frameshifts:Set[vep.VEPVariantRecord]=None):
        """ Construct a PeptideNode object.
        
        Args:
            seq (aa.AminoAcidSeqRecord): The amino acid sequence.
            variants (List[VariantRecordWithCoordinate]): The variant records
                carried in the sequence.
            frameshifts (Set[VEPVariantRecord]): Frameshifting variants.
            in_nodes (Set[PeptideNode]): Inbound nodes.
            ou_nodes (Set[PeptideNode]): Outbound nodes
        """
        self.seq = seq
        self.variants = variants
        self.in_nodes = set() if in_nodes is None else in_nodes
        self.out_nodes = set() if out_nodes is None else out_nodes
        self.frameshifts = set() if frameshifts is None else frameshifts
    
    def add_out_edge(self, node:PeptideNode) -> None:
        """ Add a outbound edge from this node.
        
        Args:
            node (PeptideNode): The outboud node of the edge to add.
        """
        self.out_nodes.add(node)
        node.in_nodes.add(self)
    
    def remove_out_edge(self, node:PeptideNode) -> None:
        """ Remove a outbound edge. It tries to remove the node, but won't
        raise any error if it doesn't exist.
        
        Args:
            node (PeptideNode): The outbound node of the edge to remove.
        """
        try:
            self.out_nodes.remove(node)
        except KeyError:
            pass
        try:
            node.in_nodes.remove(self)
        except KeyError:
            pass
        return
        
    def find_reference_next(self) -> PeptideNode:
        """ Find and return the next reference node. """
        if not self.out_nodes:
            return None
        for node in self.out_nodes:
            if not node.variants:
                return node
        raise ValueError('None of the out nodes is reference.')
    
    def find_reference_prev(self) -> PeptideNode:
        """ Find and return the previous reference nocd. """
        if not self.out_nodes:
            return None
        for node in self.in_nodes:
            if not node.variants:
                return node
        raise ValueError('None of the in nodes is reference.')
    
    def split_node(self, index:int) -> PeptideNode:
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
        for variant in self.variants:
            if variant.location.start < index:
                left_variants.append(variant)
            if variant.location.end > index:
                right_variants.append(variant.shift(-index))
        self.seq = left_seq
        self.variants = left_variants
        
        new_node = PeptideNode(seq=right_seq, variants=right_variants)
        new_node.frameshifts = copy.copy(self.frameshifts)

        while self.out_nodes:
            node = self.out_nodes.pop()
            self.remove_out_edge(node)
            new_node.add_out_edge(node)
        self.add_out_edge(new_node)
        return new_node