""" Variant Record Pool """
from __future__ import annotations
from typing import Dict, IO, Iterable, List, TYPE_CHECKING
from moPepGen import ERROR_INDEX_IN_INTRON
from . import VariantRecord, io


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.dna import DNASeqDict

T = Dict[str, List[VariantRecord]]

class VariantRecordPool():
    """ This is used to store all variant records read in from GVF files.
    Variant records (VariantRecord) are first filtered based on their location
    and put into different buckets of the pool.

    Attributes:
        transcriptional (Dict[str, List[VariantRecord]]): Variant records
            located in the transcript regions (exons), in transcriptional
            coordinates.
        intronic (Dict[str, List[VariantRecord]]): Variant records located in
            the intronic regions. In gene coordinates.
        genetic (Dict[str, List[VariantRecord]]): Variant records without a
            transcript ID (e.g., UTR). In gene coordinates.
    """
    def __init__(self, transcriptional:T=None, intronic:T=None, genetic:T=None,
            fusion:List[VariantRecord]=None):
        """ Constructor

        Args:
            transcriptional (Dict[str, List[VariantRecord]]): Variant records
                located in the transcript regions (exons), in transcriptional
                coordinates.
            intronic (Dict[str, List[VariantRecord]]): Variant records located
                in the intronic regions. In gene coordinates.
            genetic (Dict[str, List[VariantRecord]]): Variant records without a
                transcript ID (e.g., UTR). In gene coordinates.
        """
        self.transcriptional = transcriptional or {}
        self.intronic = intronic or {}
        self.genetic = genetic or {}
        self.fusion = fusion or []


    def add_genetic_variant(self, record:VariantRecord, gene_id:str=None):
        """ Add a variant with genetic coordinate """
        if not gene_id:
            gene_id = record.location.seqname

        if gene_id not in self.genetic:
            self.genetic[gene_id] = [record]
        else:
            self.genetic[gene_id].append(record)

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

    def add_fusion_variant(self, record:VariantRecord):
        """ Add a fuction variant """
        self.fusion.append(record)

    def load_variants(self, handle:IO, anno:GenomicAnnotation,
            genome:DNASeqDict):
        """ Load variants """
        for record in io.parse(handle):
            if not record.has_transcript():
                gene_id = record.location.seqname
                self.add_genetic_variant(record, gene_id)
                continue

            tx_id = record.attrs['TRANSCRIPT_ID']
            if record.is_fusion():
                record.shift_breakpoint_to_closest_exon(anno)
                tx_record = record.to_transcript_variant(anno, genome, tx_id)
                self.add_fusion_variant(tx_record)
                continue

            if record.is_spanning_over_splicing_site(anno, tx_id):
                self.add_genetic_variant(record, tx_id)
                continue

            try:
                tx_record = record.to_transcript_variant(anno, genome, tx_id)
                self.add_transcriptional_variant(tx_record, tx_id)
            # except err.FusionBreakpointIsEndOfTranscript as e:
            #     continue
            except ValueError as e:
                if e.args[0] == ERROR_INDEX_IN_INTRON:
                    self.add_intronic_variant(record, tx_id)
                else:
                    raise e

    def sort(self):
        """ sort """
        for val in self.genetic.values():
            val.sort()
        for val in self.intronic.values():
            val.sort()
        for val in self.transcriptional.values():
            val.sort()

    def filter_variants(self, anno:GenomicAnnotation, gene_id:str=None,
            tx_ids:List[str]=None, exclude_type:List[str]=None, start:int=None,
            end:int=None, intron:bool=True,
            segments:Iterable[VariantRecord]=None, return_coord:str='gene'
            ) -> List[VariantRecord]:
        """ Filter variants of located at a given position (start and end) of
        a gene.

        Arguments:
            gene_id (str): Gene ID
            start (int): Start position with genetic coordinates.
            end (int): End position with genetic coordinates. If
            exclude_type (List[str]): Variant types that should be excluded.
            anno (GenomicAnnotation): The genomic annotation object.
        """
        if return_coord not in ['gene', 'transcript']:
            raise ValueError(
                "Don't know how to return variants in coordinate of "
                f"{return_coord}. return_coord must be either 'gene' or "
                "'transcript'"
            )
        if return_coord == 'transcript' and intron is True:
            raise ValueError(
                "Don't know how to return intronic variants in transcript"
                "coordinates."
            )
        if segments:
            def _filter(x):
                for segment in segments:
                    if segment.location.is_superset(x.location):
                        return True
                return False
        elif start is not None and end is not None:
            _filter = lambda x: x.location.start > start and x.location.end < end
        elif start is not None:
            _filter = lambda x: x.location.start > start
        elif end is not None:
            _filter = lambda x: x.location.end < end
        else:
            raise ValueError('Arguments provided do not match requirement.')

        records = set()
        if tx_ids:
            _tx_ids = set(tx_ids)
        else:
            _tx_ids = set()
        if gene_id:
            _tx_ids.update(anno.genes[gene_id].transcripts)
        for tx_id in _tx_ids:
            gene_id = anno.transcripts[tx_id].transcript.gene_id
            if tx_id in self.transcriptional:
                for record in self.transcriptional[tx_id]:
                    if record.type in exclude_type:
                        continue
                    record_gene = anno.variant_coordinates_to_gene(record, gene_id)
                    if _filter(record_gene):
                        if return_coord == 'gene':
                            records.add(record_gene)
                        else:
                            records.add(record)
            if not intron:
                continue
            if tx_id in self.intronic:
                for record in self.intronic[tx_id]:
                    if record.type in exclude_type:
                        continue
                    if _filter(record):
                        if return_coord == 'gene':
                            records.add(record)
        records = list(records)
        records.sort()
        return records
