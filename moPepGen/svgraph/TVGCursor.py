""" Module for cursor for codon finding """
from moPepGen import svgraph


class TVGCursor():
    """ Cursor class, used when walking in the DNA transcirpt graph for
    finding codons. """
    def __init__(self, node:svgraph.DNANode, search_orf:bool=False):
        """ Constructor """
        self.node = node
        self.search_orf = search_orf
