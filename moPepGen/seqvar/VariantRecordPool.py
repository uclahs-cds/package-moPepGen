""" Variant Record Pool """
from __future__ import annotations
from typing import Dict, Iterable, List, TYPE_CHECKING
from pathlib import Path
from moPepGen import logger, ERROR_INDEX_IN_INTRON
from . import VariantRecord, io


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.dna import DNASeqDict

T = Dict[str, List[VariantRecord]]

class VariantRecordPool():
    """ Variant Record Pool """
    def __init__(self, transcriptional:T=None, intronic:T=None, genetic:T=None):
        """ Constructor """
        self.transcriptional = transcriptional or {}
        self.intronic = intronic or {}
        self.genetic = genetic or {}

    def add_genetic_variant(self, record:VariantRecord, gene_id:str=None):
        """ Add a variant with genetic coordinate """
        if not gene_id:
            gene_id = record.location.seqname

        if gene_id not in self.genetic:
            self.genetic[gene_id] = [record]
        else:
            self.genetic[gene_id].append(record)

    def add_genetic_variants(self, records:Iterable[VariantRecord], gene_id:str=None):
        """ Add multiple VariantRecord with genetic coordinate  """
        for record in records:
            self.add_genetic_variant(record=record, gene_id=gene_id)

    def add_intronic_variant(self, record:VariantRecord, tx_id:str):
        """ Add a variant with genetic coordinate that is in the intron of
        a transcript """
        if tx_id not in self.intronic:
            self.intronic[tx_id] = [record]
        else:
            self.intronic[tx_id].append(record)

    def add_intronic_variants(self, records:Iterable[VariantRecord], tx_id:str):
        """ Add multiple intronic VariantRecord with genetic coordinate """
        for record in records:
            self.add_intronic_variant(record=record, tx_id=tx_id)

    def add_transcriptional_variant(self, record:VariantRecord, tx_id:str=None):
        """ Add a variant with transcriptional coordinate """
        if not tx_id:
            tx_id = record.location.seqname

        if tx_id not in self.transcriptional:
            self.transcriptional[tx_id] = [record]
        else:
            self.transcriptional[tx_id].append(record)

    def add_transcriptional_variants(self, records:Iterable[VariantRecord],
            tx_id:str=None):
        """ Add multiple VariantRecords with transcriptional coordinates """
        for record in records:
            self.add_transcriptional_variant(record, tx_id)

    @staticmethod
    def load_variants(files:Path, anno:GenomicAnnotation,
            genome:DNASeqDict, verbose:bool) -> VariantRecordPool:
        """ Load variants from files """
        variants = VariantRecordPool()
        for file in files:
            for record in io.parse(file):
                if not record.has_transcripts():
                    gene_id = record.location.seqname
                    variants.add_genetic_variant(record, gene_id)
                    continue

                for tx_id in record.attrs['TRANSCRIPTS']:
                    try:
                        tx_records = record.to_transcript_variants(anno, genome, [tx_id])
                    except ValueError as e:
                        if e.args[0] == ERROR_INDEX_IN_INTRON:
                            variants.add_intronic_variant(record, tx_id)
                            continue
                    variants.add_transcriptional_variants(tx_records, tx_id)
            if verbose:
                logger(f'Variant file {file} loaded.')

        for val in variants.genetic.values():
            val.sort()
        for val in variants.intronic.values():
            val.sort()
        for val in variants.transcriptional.values():
            val.sort()

        if verbose:
            logger('Variant records sorted.')

        return variants
