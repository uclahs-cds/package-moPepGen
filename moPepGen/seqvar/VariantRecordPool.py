""" Variant Record Pool """
from __future__ import annotations
import copy
from typing import Dict, IO, Iterable, List, TYPE_CHECKING, Union
from moPepGen import ERROR_INDEX_IN_INTRON, circ
from . import VariantRecord, io, GVFMetadata
from .VariantRecordPoolOnDisk import TranscriptionalVariantSeries


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
    def __init__(self, data:Dict[str, TranscriptionalVariantSeries]=None,
            anno:GenomicAnnotation=None):
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
        self.data = data or {}
        self.anno = anno

    def __contins__(self, key:str):
        """ in """
        return key in self.data

    def __getitem__(self, key:str) -> Union[TranscriptionalVariantSeries]:
        """ get item """
        return self.data[key]

    def __setitem__(self, key:str, val:Union[TranscriptionalVariantSeries]):
        """ setitem """
        self.data[key] = val

    def __iter__(self) -> Iterable[str]:
        """ generator """
        for key in self.data:
            yield key

    def __copy__(self) -> VariantRecordPool:
        """ copy """
        data = copy.copy(self.data)
        anno = copy.copy(self.anno)
        return self.__class__(data, anno)

    def add_intronic_variant(self, record:VariantRecord, tx_id:str=None):
        """ Add a variant with genetic coordinate that is in the intron of
        a transcript """
        if tx_id is None:
            tx_id = record.transcript_id
        if tx_id not in self:
            self[tx_id] = TranscriptionalVariantSeries()
        self[tx_id].intronic.append(record)

    def add_intronic_variants(self, records:Iterable[VariantRecord], tx_id:str):
        """ Add multiple intronic VariantRecord with genetic coordinate """
        for record in records:
            self.add_intronic_variant(record=record, tx_id=tx_id)

    def add_transcriptional_variant(self, record:VariantRecord, tx_id:str=None):
        """ Add a variant with transcriptional coordinate """
        if tx_id is None:
            tx_id = record.location.seqname
        if tx_id not in self:
            self[tx_id] = TranscriptionalVariantSeries()
        self[tx_id].transcriptional.append(record)

    def add_fusion_variant(self, record:VariantRecord, tx_id:str=None):
        """ Add a fuction variant """
        if tx_id is None:
            tx_id = record.location.seqname
        if tx_id not in self:
            self[tx_id] = TranscriptionalVariantSeries()
        self[tx_id].fusion.append(record)

    def add_circ_rna(self, record:circ.CircRNAModel, tx_id:str=None):
        """ Add a circRNA """
        if tx_id is None:
            tx_id = record.transcript_id
        if tx_id not in self:
            self[tx_id] = TranscriptionalVariantSeries()
        self[tx_id].circ_rna.append(record)

    def load_variants(self, handle:IO, anno:GenomicAnnotation,
            genome:DNASeqDict):
        """ Load variants """
        handle.seek(0)
        metadata = GVFMetadata.parse(handle)
        handle.seek(0)

        if metadata.is_circ_rna():
            for record in circ.io.parse(handle):
                self.add_circ_rna(record)
            return

        for record in io.parse(handle):
            tx_id = record.transcript_id
            if record.is_fusion():
                record.shift_breakpoint_to_closest_exon(anno)
                tx_record = record.to_transcript_variant(anno, genome, tx_id)
                self.add_fusion_variant(tx_record, tx_id)
                continue

            if record.is_spanning_over_splicing_site(anno, tx_id):
                continue

            try:
                tx_record = record.to_transcript_variant(anno, genome, tx_id)
                if tx_record.type == 'Deletion':
                    tx_model = anno.transcripts[tx_id]
                    chrom = tx_model.transcript.chrom
                    tx_seq = tx_model.get_transcript_sequence(genome[chrom])
                    tx_record.shift_deletion_up(tx_seq)
                self.add_transcriptional_variant(tx_record, tx_id)

            except ValueError as e:
                if e.args[0] == ERROR_INDEX_IN_INTRON:
                    self.add_intronic_variant(record, tx_id)
                else:
                    raise e
        self.sort()

    def sort(self):
        """ sort """
        for key in self:
            series = self[key]
            if isinstance(series, TranscriptionalVariantSeries):
                series.sort()

    def filter_variants(self, gene_id:str=None, tx_ids:List[str]=None,
            exclude_type:List[str]=None, start:int=None, end:int=None,
            intron:bool=True, segments:Iterable[VariantRecord]=None,
            return_coord:str='gene') -> List[VariantRecord]:
        """ Filter variants of located at a given position (start and end) of
        a gene.

        Arguments:
            gene_id (str): Gene ID
            start (int): Start position with genetic coordinates.
            end (int): End position with genetic coordinates. If
            exclude_type (List[str]): Variant types that should be excluded.
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
                    if segment.location.start < x.location.start < \
                            x.location.end < segment.location.end:
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
            _tx_ids.update(self.anno.genes[gene_id].transcripts)
        for tx_id in _tx_ids:
            gene_id = self.anno.transcripts[tx_id].transcript.gene_id
            if tx_id in self:
                for record in self[tx_id].transcriptional:
                    if record.type in exclude_type:
                        continue
                    try:
                        record_gene = self.anno.variant_coordinates_to_gene(record, gene_id)
                    except ValueError as e:
                        if record.is_merged_mnv():
                            continue
                        raise e
                    if _filter(record_gene):
                        if return_coord == 'gene':
                            records.add(record_gene)
                        else:
                            records.add(record)
            if not intron:
                continue
            if tx_id in self:
                for record in self[tx_id].intronic:
                    if record.type in exclude_type:
                        continue
                    if _filter(record):
                        if return_coord == 'gene':
                            records.add(record)
        records = list(records)
        records.sort()
        return records
