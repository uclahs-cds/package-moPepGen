""" Module for cursor for codon finding """
from __future__ import annotations
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from moPepGen.svgraph.TVGNode import TVGNode

class TVGCursor():
    """ Cursor class, used when walking in the DNA transcirpt graph for
    finding codons. """
    def __init__(self, node:TVGNode, search_orf:bool=False):
        """ Constructor """
        self.node = node
        self.search_orf = search_orf
