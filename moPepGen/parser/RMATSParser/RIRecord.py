""" Module for Retained Intron """
from .RMATSRecord import RMATSRecord


class RIRecord(RMATSRecord):
    """ Retained Intron """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            first_exon_start:int, first_exon_end:int, second_exon_start:int,
            second_exon_end:int, upstream_exon_start:int,
            upstream_exon_end:int, downstream_exon_start:int,
            downstream_exon_end:int):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.first_exon_start = first_exon_start
        self.first_exon_end = first_exon_end
        self.second_exon_start = second_exon_start
        self.second_exon_end = second_exon_end
        self.upstream_exon_start = upstream_exon_start
        self.upstream_exon_end = upstream_exon_end
        self.downstream_exon_start = downstream_exon_start
        self.downstream_exon_end = downstream_exon_end
