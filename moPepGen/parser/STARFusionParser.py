""" Module for STAR-Fusion parser """
from __future__ import annotations
from typing import Iterable, List, Dict
import itertools
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, gtf, dna


GeneAnnotations = Dict[str, Dict[str, gtf.TranscriptAnnotationModel]]

def parse(path:str) -> Iterable[STARFusionRecord]:
    """ Parse the STAR-Fusion's output and returns an iterator.

    Args:
        path (str): Path to the STAR-Fusion's output file.
    """
    with open(path, 'r') as handle:
        line = next(handle, None)
        while line:
            if line.startswith('#'):
                line = next(handle, None)
                continue
            fields = line.rstrip().split('\t')
            yield STARFusionRecord(
                fusion_name=fields[0],
                junction_read_count=int(fields[1]),
                spanning_frag_count=int(fields[2]),
                est_j=float(fields[3]),
                est_s=float(fields[4]),
                splice_type=fields[5],
                left_gene=fields[6].split('^')[-1],
                left_breakpoint=fields[7],
                right_gene=fields[8].split('^')[-1],
                right_breakpoint=fields[9],
                junction_reads=fields[10].split(','),
                spanning_frags=fields[11].split(','),
                large_anchor_support=fields[12],
                ffpm=float(fields[13]),
                left_break_dinuc=fields[14],
                left_break_entropy=float(fields[15]),
                right_break_dinuc=fields[16],
                right_break_entropy=float(fields[17]),
                annots=fields[18].strip('][').replace('"', '').split(',')
            )
            line = next(handle, None)


class STARFusionRecord():
    """ Defines a record of STAR-Fusion's output. """
    def __init__(self, fusion_name:str, junction_read_count:int,
            spanning_frag_count:int, est_j:float, est_s:float,
            splice_type:str, left_gene:str, left_breakpoint:str,
            right_gene:str, right_breakpoint:str, junction_reads:List[str],
            spanning_frags:List[str], large_anchor_support:str, ffpm:float,
            left_break_dinuc:str, left_break_entropy:float,
            right_break_dinuc:str, right_break_entropy:float, annots:List[str]
            ):
        """"""
        self.fusion_name = fusion_name
        self.junction_read_count = junction_read_count
        self.spanning_frag_count = spanning_frag_count
        self.est_j = est_j
        self.est_s = est_s
        self.splice_type = splice_type
        self.left_gene = left_gene
        self.left_breakpoint = left_breakpoint
        self.right_gene = right_gene
        self.right_breakpoint = right_breakpoint
        self.junction_reads = junction_reads
        self.spanning_frags = spanning_frags
        self.large_anchor_support = large_anchor_support
        self.ffpm = ffpm
        self.left_break_dinuc = left_break_dinuc
        self.left_break_entropy = left_break_entropy
        self.right_break_dinuc = right_break_dinuc
        self.right_break_entropy = right_break_entropy
        self.annots = annots

    def convert_to_variant_records(self, anno:GeneAnnotations,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        """ COnvert a STAR-Fusion's record to VariantRecord.

        Args:
            anno (Dict[str, Dict[str, gtf.TranscriptAnnotationModel]]): The
                genome annotation. The first level key is the gene ID, and the
                second level key is the transcript ID.
            genome (dna.DNASeqDict): The genome sequence records.
        Returns:
            List of VariantRecord
        """
        donor_transcripts = anno[self.left_gene]
        accepter_transcripts = anno[self.right_gene]
        records = []

        perms = itertools.product(donor_transcripts.keys(), accepter_transcripts.keys())
        for donor_id, accepter_id in perms:
            donor_model = donor_transcripts[donor_id]
            donor_gene_symbol = donor_model.transcript.attributes['gene_name']
            left_breakpoint = int(self.left_breakpoint.split(':')[1])
            donor_chrom = self.left_breakpoint.split(':')[0]
            donor_position = donor_model.get_transcript_index(left_breakpoint)
            seq = donor_model.get_transcript_sequence(genome[donor_chrom])
            donor_genome_position = f'{donor_chrom}:{left_breakpoint}:{left_breakpoint}'

            accepter_model = accepter_transcripts[accepter_id]
            right_breakpoint = int(self.right_breakpoint.split(':')[1])
            accepter_position = accepter_model.get_transcript_index(right_breakpoint)
            accepter_gene_symbol = accepter_model.transcript.attributes['gene_name']
            accepter_chrom = self.right_breakpoint.split(':')[0]
            accepter_genome_position = f'{accepter_chrom}:{right_breakpoint}:{right_breakpoint}'

            location = FeatureLocation(
                seqname=donor_id,
                start=donor_position,
                end=donor_position + 1
            )
            _id = f'FUSION_{donor_id}:{donor_position}-{accepter_id}:{accepter_position}'
            attrs = {
                'GENE_ID': self.left_gene,
                'GENE_SYMBOL': donor_gene_symbol,
                'GENOMIC_POSITION': donor_genome_position,
                'ACCEPTER_GENE_ID': self.right_gene,
                'ACCEPTER_TRANSCRIPT_ID': accepter_id,
                'ACCEPTER_SYMBOL': accepter_gene_symbol,
                'ACCEPTER_POSITION': accepter_position,
                'ACCEPTER_GENOMIC_POSITION': accepter_genome_position
            }
            record = seqvar.VariantRecord(
                location=location,
                ref=seq[donor_position],
                alt='<FUSION>',
                _type='Fusion',
                _id=_id,
                attrs=attrs
            )
            records.append(record)
        return records
