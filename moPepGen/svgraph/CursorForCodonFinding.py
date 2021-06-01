""""""
from moPepGen import svgraph


class CursorForCodonFinding():
    """"""
    def __init__(self, node:svgraph.DNANode, position_to_orf:str='up',
            search_orf:bool=False, start_codon_index:int=None,
            stop_codon_index:int=None):
        """"""
        self.node = node
        self.position_to_orf = position_to_orf
        self.search_orf = search_orf
        self.start_codon_index = start_codon_index
        self.stop_codon_index = stop_codon_index