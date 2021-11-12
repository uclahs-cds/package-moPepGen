""" Module for rMATS parser """
from __future__ import annotations
from typing import List
from abc import ABC, abstractmethod
from moPepGen import seqvar

class RMATSRecord(ABC):
    """ rMATS record """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str):
        """"""
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol
        self.chrom = chrom

    @abstractmethod
    def convert_to_variant_records(self, anno, genome, min_ijc, min_sjc
            ) -> List[seqvar.VariantRecord]:
        """ Convert to variant records """
