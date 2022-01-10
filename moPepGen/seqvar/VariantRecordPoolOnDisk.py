""" Variant Record Pool """
from __future__ import annotations
from typing import Dict, IO, Iterable, List, TYPE_CHECKING, Union
from pathlib import Path
from moPepGen import ERROR_INDEX_IN_INTRON, check_sha512
from moPepGen.circ.CircRNA import CircRNAModel
from moPepGen.seqvar.GVFIndex import GVFPointer, iterate_pointer
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from . import VariantRecord


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.dna import DNASeqDict


class TranscriptionalVariantSeries():
    """ Variants associated with a particular transcript """
    def __init__(self, transcriptional:List[VariantRecord]=None,
            intronic:List[VariantRecord]=None, fusion:List[VariantRecord]=None,
            circ_rna:List[CircRNAModel]=None
            ) -> None:
        """ constructor """
        self.transcriptional = transcriptional or []
        self.intronic = intronic or []
        self.fusion = fusion or []
        self.circ_rna = circ_rna or []

    def sort(self):
        """ sort """
        self.transcriptional.sort()
        self.intronic.sort()
        self.fusion.sort()

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

    def __enter__(self):
        """ enter """
        for file in self.gvf_files:
            gvf_handle = file.open('rb')
            self.gvf_handles.append(gvf_handle)
            idx_path = file.with_suffix(file.suffix + '.idx')
            if idx_path.exists():
                self.validate_gvf_index(file, idx_path)
                self.load_index(idx_path, file, gvf_handle)
            else:
                self.generate_index(file, gvf_handle)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        """ """
        for handle in self.gvf_handles:
            handle.close()

    def __contains__(self, key:str):
        """ conteins """
        return key in self.pointers

    def __iter__(self) -> Iterable[str]:
        """ generator """
        for key in self.pointers:
            yield key

    def __getitem__(self, key:str) -> TranscriptionalVariantSeries:
        """ Load variants and return as a TranscriptVariants object """
        records:List[Union[VariantRecord, CircRNAModel]] = []
        for pointer in self.pointers[key]:
            records += pointer.load()
        series = TranscriptionalVariantSeries()
        for record in records:
            if isinstance(record, CircRNAModel):
                series.circ_rna.append(record)
                continue

            tx_id = record.transcript_id

            if record.is_fusion():
                record.shift_breakpoint_to_closest_exon(self.anno)
                tx_record = record.to_transcript_variant(
                    self.anno, self.genome, tx_id
                )
                series.fusion.append(tx_record)
                continue

            if record.is_spanning_over_splicing_site(self.anno, tx_id):
                continue

            try:
                tx_record = record.to_transcript_variant(
                    self.anno, self.genome, tx_id
                )
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
