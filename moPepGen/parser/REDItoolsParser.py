""" Module for REDItoolsParser """
from __future__ import annotations
from typing import List, Tuple, Iterable
import re
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen import seqvar, gtf, ERROR_INDEX_IN_INTRON


def parse(path:str, transcript_id_column:int=16
        ) -> Iterable[REDItoolsRecord]:
    """ Parse the REDItools output table and returns an iterator. The input
    table must be annotated using the AnnotateTable.py script from REDItools 1.

    Args:
        path (str): Path to the REDItools output table.
        transcript_id_column (int): The column that holds the transcript_id.
            Defaults to 10.

    Return:
        A iterable of REDItoolsRecord.
    """
    with open(path, 'r') as handle:
        line = next(handle, None)
        line = next(handle, None)
        while line:
            line = line.rstrip()
            line = re.sub('[,&$]$', '', line)
            fields = line.split('\t')
            base_count = [int(i) for i in fields[6].strip('][')\
                    .split(', ')]
            all_subs = []
            for sub in fields[7].split(' '):
                if len(sub) > 2:
                    raise ValueError('length of sub larger than 2')
                all_subs.append((sub[0], sub[1]))

            try:
                g_coverage_q = int(fields[9])
            except ValueError:
                g_coverage_q = None

            if len(fields) <= transcript_id_column:
                raise ValueError('transcript_id_column invalid.')
            _data = re.split(r',|&|\$', fields[transcript_id_column])
            transcript_ids = [tuple(x.split('-')) for x in _data]

            yield REDItoolsRecord(
                region=fields[0],
                position=int(fields[1]),
                reference=fields[2],
                strand=int(fields[3]),
                coverage_q=int(fields[4]),
                mean_quality=float(fields[5]),
                base_count=base_count,
                all_subs=all_subs,
                frequency=float(fields[8]),
                g_coverage_q=g_coverage_q,
                transcript_id=transcript_ids
            )
            line = next(handle, None)


class REDItoolsRecord():
    """ A REDItoolsRecord defines the attributes from a REDItools output table.

    Attributes:
        region (str): chromosome sequence name.
        position (int): 1 based position.
        reference (str): nucleotide in reference.
        strand (int): 1 for position; -1 for negative.
        coverage_q30 (int): depth per site at quality score of 30.
        mean_quality (int): mean quality score per site.
        base_cout (List[int]): base distribution per site in the order of
            A, C, G and T
        all_subs (List[Tuple[str, str]]): all substritutions.
        frequency (int): observed frequency of substitution.
        transcirpt_id (List[Tuple[str,str]]): transcript IDs.
    """
    def __init__(self, region:str, position:int, reference:str, strand:int,
            coverage_q:int, mean_quality:float, base_count:List[int],
            all_subs:List[Tuple[str, str]], frequency:int, g_coverage_q:int,
            transcript_id:List[Tuple[str, str]]):
        """ constructor """
        self.region = region
        self.position = position
        self.reference = reference
        self.strand = strand
        self.all_subs = all_subs
        self.coverage_q = coverage_q
        self.mean_quality = mean_quality
        self.base_count = base_count
        self.all_subs = all_subs
        self.frequency = frequency
        self.g_coverage_q = g_coverage_q
        self.transcript_id = transcript_id
        self.base_count_order = {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3
        }

    def get_valid_subs(self, min_coverage_alt:int, min_frequency_alt:float,
            min_coverage_rna:int, min_coverage_dna:int) -> bool:
        """ Get all valid substitutions. """
        total_count = sum(self.base_count)
        valid_subs = []

        if total_count < min_coverage_rna:
            return valid_subs

        if self.g_coverage_q != -1:
            if self.g_coverage_q is None or self.g_coverage_q < min_coverage_dna:
                return valid_subs

        for sub in self.all_subs:
            alt = sub[1]
            read_count = self.base_count[self.base_count_order[alt]]
            if read_count < min_coverage_alt:
                continue
            if read_count / total_count < min_frequency_alt:
                continue
            valid_subs.append(sub)
        return valid_subs

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            min_coverage_alt:int, min_frequency_alt:float,
            min_coverage_rna:int, min_coverage_dna:int
            ) -> List[seqvar.VariantRecord]:
        """ Convert to VariantRecord.

        Args:
            anno (gtf.TranscriptFTFDict): the reference genome annotations.

        Return:
            A list of seqvar.VariantRecord. Multiple variant records can be
            returned, because a genomic position can be on multiple transcripts
            with different splicing.
        """
        _ids = []
        for tx_id, feature in self.transcript_id:
            if feature == 'transcript':
                _ids.append(tx_id)

        records = []
        genomic_location = f"{self.region}:{self.position}"
        for tx_id in _ids:
            tx_model:gtf.TranscriptAnnotationModel = anno.transcripts[tx_id]
            try:
                position = tx_model.get_transcript_index(self.position - 1)
            except ValueError as e:
                if e.args[0] == ERROR_INDEX_IN_INTRON:
                    continue
            gene_id = tx_model.transcript.gene_id
            gene_model = anno.genes[gene_id]
            strand = gene_model.strand
            position = anno.coordinate_genomic_to_gene(self.position - 1, gene_id)
            location = FeatureLocation(
                seqname=gene_id,
                start=position,
                end=position + 1
            )
            valid_subs = self.get_valid_subs(
                min_coverage_alt=min_coverage_alt,
                min_frequency_alt=min_frequency_alt,
                min_coverage_rna=min_coverage_rna,
                min_coverage_dna=min_coverage_dna
            )
            for sub in valid_subs:
                ref = sub[0]
                alt = sub[1]

                _id = f'RES-{position + 1}-{ref}-{alt}'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'GENOMIC_POSITION': genomic_location,
                    'STRAND': strand
                }
                record = VariantRecord(
                    location=location,
                    ref=ref,
                    alt=alt,
                    _type='RNAEditingSite',
                    _id=_id,
                    attrs=attrs
                )
                records.append(record)
        return records
