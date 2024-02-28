""" Variant Record Pool """
from __future__ import annotations
import copy
from typing import Dict, IO, Iterable, List, TYPE_CHECKING, Union
from pathlib import Path
from moPepGen import ERROR_INDEX_IN_INTRON, check_sha512, circ, constant
from moPepGen.seqvar.GVFIndex import GVFPointer, iterate_pointer
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from . import VariantRecord


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.dna import DNASeqDict, DNASeqRecordWithCoordinates


class TranscriptionalVariantSeries():
    """ Variants associated with a particular transcript """
    def __init__(self, transcriptional:List[VariantRecord]=None,
            intronic:List[VariantRecord]=None, fusion:List[VariantRecord]=None,
            circ_rna:List[circ.CircRNAModel]=None
            ) -> None:
        """ constructor """
        self.transcriptional = transcriptional or []
        self.intronic = intronic or []
        self.fusion = fusion or []
        self.circ_rna = circ_rna or []

    def __copy__(self) -> TranscriptionalVariantSeries:
        """ copy """
        return self.__class__(
            transcriptional=copy.copy(self.transcriptional),
            intronic=copy.copy(self.intronic),
            fusion=copy.copy(self.fusion),
            circ_rna=copy.copy(self.circ_rna)
        )

    def sort(self):
        """ Sort each slot of variants in order """
        self.transcriptional.sort()
        self.intronic.sort()
        self.fusion.sort()

    def get_additional_transcripts(self) -> List[str]:
        """ Get any additioanl transcript IDs that are associated with any
        variants, for example, fusion. """
        transcripts = []
        for variant in self.fusion:
            transcripts.append(variant.attrs['ACCEPTER_TRANSCRIPT_ID'])
        return transcripts

    def is_gene_sequence_needed(self) -> bool:
        """ Check if the gene sequence is needed when calling variant peptides.
        Usually for variants that invole retaining of a complete or partial
        intron. """
        return len(self.fusion) > 0 or len(self.circ_rna) > 0 or \
            any(x.type in ['Insertion', 'Substitution'] for x in self.transcriptional)

    def is_empty(self) -> bool:
        """ check if the series is empty """
        return len(self.transcriptional) == 0 and len(self.fusion) == 0 and \
            len(self.circ_rna) == 0

    def has_any_noncanonical_transcripts(self) -> bool:
        """ check if the series has any noncanonical transcripts """
        return len(self.fusion) > 0 or len(self.circ_rna) > 0 or \
            any(x.type in constant.ALTERNATIVE_SPLICING_TYPES for x in self.transcriptional)

    def has_any_alternative_splicing(self) -> bool:
        """ Check if there is any alternative splicing """
        return any(x.type in constant.ALTERNATIVE_SPLICING_TYPES for x in self.transcriptional)

    def get_highest_hypermutated_region_complexity(self, distance:int=6):
        """ Calculate the number of variants in the most hypermutated region """
        end = 0
        max_n = 0
        cur_n = 0
        for variant in self.transcriptional:
            if variant.type in constant.ALTERNATIVE_SPLICING_TYPES:
                continue
            if variant.location.start - end >= distance:
                end = variant.location.end
                max_n = max(max_n, cur_n)
                cur_n = 1
                continue
            cur_n += 1
            end = variant.location.end
        return max(max_n, cur_n)

class VariantRecordPoolOnDisk():
    """ Variant record pool in disk """
    def __init__(self, pointers:Dict[str, List[GVFPointer]]=None,
            gvf_files:List[Path]=None, gvf_handles:List[IO]=None,
            anno:GenomicAnnotation=None, genome:DNASeqDict=None) -> None:
        """ constructor """
        self.pointers = pointers or {}
        self.gvf_files = gvf_files or []
        self.gvf_handles = gvf_handles or []
        self.anno = anno
        self.genome = genome

    def __contains__(self, key:str):
        """ conteins """
        return key in self.pointers

    def __iter__(self) -> Iterable[str]:
        """ generator """
        for key in self.pointers:
            yield key

    def __getitem__(self, key:str) -> TranscriptionalVariantSeries:
        """ Load variants and return as a TranscriptVariants object """
        records:List[Union[VariantRecord, circ.CircRNAModel]] = []
        for pointer in self.pointers[key]:
            records += pointer.load()
        records = set(records)
        series = TranscriptionalVariantSeries()
        cached_seqs:Dict[str, DNASeqRecordWithCoordinates] = {}
        for record in records:
            if isinstance(record, circ.CircRNAModel):
                series.circ_rna.append(record)
                continue

            tx_id = record.transcript_id

            if record.is_fusion():
                record.shift_breakpoint_to_closest_exon(self.anno)
                tx_record = record.to_transcript_variant(
                    self.anno, self.genome, tx_id, cached_seqs
                )
                series.fusion.append(tx_record)
                continue

            if record.is_spanning_over_splicing_site(self.anno, tx_id):
                continue

            try:
                tx_record = record.to_transcript_variant(
                    self.anno, self.genome, tx_id
                )
                if tx_record.type == 'Deletion':
                    if tx_id in cached_seqs:
                        tx_seq = cached_seqs[tx_id]
                    else:
                        tx_model = self.anno.transcripts[tx_id]
                        chrom = tx_model.transcript.chrom
                        tx_seq = tx_model.get_transcript_sequence(self.genome[chrom])
                    tx_record.shift_deletion_up(tx_seq)
                series.transcriptional.append(tx_record)
            # except err.FusionBreakpointIsEndOfTranscript as e:
            #     continue
            except ValueError as e:
                if e.args[0] == ERROR_INDEX_IN_INTRON:
                    series.intronic.append(record)
                else:
                    raise e

        series.sort()
        return series

    def load_index(self, index_file:Path, gvf_file:Path, gvf_handle:IO):
        """ Load index file """
        with open(gvf_file, 'rt') as handle:
            metadata = GVFMetadata.parse(handle)
        is_circ_rna = metadata.is_circ_rna()
        gvf_handle.seek(0)
        with open(index_file, 'rt') as idx_handle:
            reader = GVFPointer.parse(
                index_handle=idx_handle, gvf_handle=gvf_handle,
                is_circ_rna=is_circ_rna
            )
            for pointer in reader:
                if pointer.key in self.pointers:
                    self.pointers[pointer.key].append(pointer)
                else:
                    self.pointers[pointer.key] = [pointer]

    def generate_index(self, gvf_file:Path, gvf_handle:IO):
        """ generate index """
        with open(gvf_file, 'rt') as handle:
            metadata = GVFMetadata.parse(handle)
        is_circ_rna = metadata.is_circ_rna()
        gvf_handle.seek(0)
        for pointer in iterate_pointer(gvf_handle, is_circ_rna):
            if pointer.key in self.pointers:
                self.pointers[pointer.key].append(pointer)
            else:
                self.pointers[pointer.key] = [pointer]
        gvf_handle.seek(0)

    @staticmethod
    def validate_gvf_index(gvf_file:Path, idx_file:Path):
        """ Validate the GVF file with index """
        sum_expect = None
        with open(gvf_file, 'rb') as gvf_handle:
            sum_actual = check_sha512(gvf_handle)

        with open(idx_file, 'rt') as idx_handle:
            for line in idx_handle:
                line:str
                if line.startswith('#'):
                    line = line.rstrip().lstrip('# ')
                    if line.startswith('CHECKSUM='):
                        sum_expect = line.split('=')[1]

                        break
                else:
                    break

        if sum_expect is None:
            raise ValueError('Cannot find checksum value from the idx file.')
        if sum_actual == sum_expect:
            return True

        raise ValueError("GVF checksum don't match.")

    def filter_variants(self, gene_id:str=None, tx_ids:List[str]=None,
            exclude_type:List[str]=None, start:int=None, end:int=None,
            intron:bool=True, segments:Iterable[VariantRecord]=None,
            return_coord:str='gene') -> List[VariantRecord]:
        """ filter variants """
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
            _tx_ids.update(self.anno.genes[gene_id].transcripts)

        for tx_id in _tx_ids:
            gene_id = self.anno.transcripts[tx_id].transcript.gene_id
            if tx_id not in self.pointers:
                continue
            series = self[tx_id]
            for record in series.transcriptional:
                if record.type in exclude_type:
                    continue
                record_gene = self.anno.variant_coordinates_to_gene(record, gene_id)
                if _filter(record_gene):
                    if return_coord == 'gene':
                        records.add(record_gene)
                    else:
                        records.add(record)

            if not intron:
                continue

            for record in series.intronic:
                if record.type in exclude_type:
                    continue
                if _filter(record):
                    if return_coord == 'gene':
                        records.add(record)

        records = list(records)
        records.sort()
        return records

    def get_transcript_order(self) -> List[str]:
        """ Get the order of transcripts to call variant peptides """
        def sort_key(key:str):
            return self[key].get_highest_hypermutated_region_complexity()
        return sorted(self.pointers.keys(), key=sort_key, reverse=True)

class VariantRecordPoolOnDiskOpener():
    """ Helper class to open all GVF files of a VariantRecordPoolOnDisk. """
    def __init__(self, pool:VariantRecordPoolOnDisk):
        """ Constructor """
        self.pool = pool

    def __enter__(self):
        """ enter """
        self.open()
        return self.pool

    def __exit__(self, exception_type, exception_value, exception_traceback):
        """ exit """
        self.close()

    def open(self):
        """ Open all GVF files """
        for file in self.pool.gvf_files:
            gvf_handle = file.open('rb')
            self.pool.gvf_handles.append(gvf_handle)
            idx_path = file.with_suffix(file.suffix + '.idx')
            if idx_path.exists():
                self.pool.validate_gvf_index(file, idx_path)
                self.pool.load_index(idx_path, file, gvf_handle)
            else:
                self.pool.generate_index(file, gvf_handle)

    def close(self):
        """ Close all file handles """
        for handle in self.pool.gvf_handles:
            handle.close()
