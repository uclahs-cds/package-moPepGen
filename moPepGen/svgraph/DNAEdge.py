""" Module for ENAEdge class """
from moPepGen import svgraph


class DNAEdge():
    """ Defines the edges in the TranscriptVariantGraph

    Attributes:
        in_node (Node): The inbond node.
        out_node (Node): The outbond node.
        type (str): The edge type. Must be either of orf_start, orf_end,
            variant_start, variant_end, cleave, or reference
    """
    def __init__(self, in_node:svgraph.DNANode, out_node:svgraph.DNANode,
            _type:str):
        """ Constructor for Edge

        Args:
            in_node (Node): The inbond node.
            out_node (Node): The outbond node.
            type (str): The edge type. Must be either of orf_start, orf_end,
                mutation_start, mutation_end, cleave, or reference
        """
        self.in_node = in_node
        self.out_node = out_node
        edge_types = ['variant_start', 'variant_end', 'reference']
        if _type not in edge_types:
            raise ValueError(f'type {_type} not from {edge_types}')
        self.type = _type
