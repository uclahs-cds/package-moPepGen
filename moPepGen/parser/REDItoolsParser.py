""" Module for REDIToolsParser """
from __future__ import annotations
from typing import List, Tuple, Iterable
import re
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen import seqvar, gtf


def parse(path:str, transcript_id_column:int=16
        ) -> Iterable[REDIToolsRecord]:
    """ Parse REDITools output table and returns an iterator. The input
    table must be annotated using the AnnotateTable.py script from REDITools 1.

    Args:
        path (str): Path to the REDITools output table.
        transcript_id_column (int): The column that holds the transcript_id.
            Defaults to 10.

    Return:
        A iterable of REDIToolsRecord.
    """
    with open(path, 'r') as handle:
        next(handle)
        for line in handle:
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

            if len(fields) <= transcript_id_column:
                raise ValueError('transcript_id_column invalid.')
            _data = re.split(',|&|\$', fields[transcript_id_column])
            transcript_ids = [tuple(x.split('-')) for x in _data]

            yield REDIToolsRecord(
                region=fields[0],
                position=int(fields[1]),
                reference=fields[2],
                strand=int(fields[3]),
                coverage_q30=int(fields[4]),
                mean_quality=float(fields[5]),
                base_count=base_count,
                all_subs=all_subs,
                frequency=float(fields[8]),
                transcript_id=transcript_ids
            )

class REDIToolsRecord():
    """ A REDIToolsRecord defines the attributes from a REDITools output table.

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
            coverage_q30:int, mean_quality:float, base_count:List[int],
            all_subs:List[Tuple[str, str]], frequency=int,
            transcript_id=List[Tuple[str, str]]):
        """ constructor """
        self.region = region
        self.position = position
        self.reference = reference
        self.starnd = strand
        self.all_subs = all_subs
        self.coverage_q30 = coverage_q30
        self.mean_quality = mean_quality
        self.base_count = base_count
        self.all_subs = all_subs
        self.frequency = frequency
        self.transcript_id = transcript_id


    def convert_to_variant_records(self, anno:gtf.TranscriptGTFDict
            ) -> List[seqvar.VariantRecord]:
        """ Convert to variant_records.

        Args:
            anno (gtf.TranscriptFTFDict): the reference genome annotations.

        Return:
            A list of seqvar.VariantRecord. Multiple variant records can be
            returned, because a genomic position can be on multiple transcripts
            with different splicing.
        """
        _ids = []
        for transcript_id, feature in self.transcript_id:
            if feature == 'transcript':
                _ids.append(transcript_id)

        records = []
        for transcript_id in _ids:
            model:gtf.TranscriptAnnotationModel = anno[transcript_id]
            try:
                position = model.get_transcript_index(self.position - 1)
            except ValueError as e:
                if e.args[0] == gtf.INDEX_IN_INTRON_ERROR:
                    continue
            location=FeatureLocation(
                seqname=transcript_id,
                start=position,
                end=position + 1
            )
            for sub in self.all_subs:
                ref = sub[0]
                alt = sub[1]
                if model.transcript.strand == -1:
                    ref = str(Seq(ref).complement())
                    alt = str(Seq(alt).complement())
                _id = f'RNA_editing_site-{ref}-{alt}'
                record = VariantRecord(
                    location=location,
                    ref=ref,
                    alt=alt,
                    _type='RNAEditingSite',
                    _id=_id
                )
                records.append(record)
        return records
