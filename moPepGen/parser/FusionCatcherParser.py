""" Module for FusionCatcher parser """
from __future__ import annotations
from typing import Iterable, List, Dict
import itertools
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, gtf, dna


GeneAnnotations = Dict[str, Dict[str, gtf.TranscriptAnnotationModel]]

def parse(path:str) -> Iterable[FusionCatcherRecord]:
    """ Parse the FusionCatcher's output and returns an iterator.

    Args:
        path (str): Path to the FusionCatcher's output file.
    """
    with open(path, 'r') as handle:
        line = next(handle, None)
        while line:
            if line.startswith('#'):
                line = next(handle, None)
                continue
            fields = line.rstrip().split('\t')
            yield FusionCatcherRecord(
                five_end_gene_symbol = fields[0],
                three_end_gene_symbol = fields[1],
                fusion_descriptions = fields[2].split(','),
                multimapping_read_count = int(fields[3]),
                junction_pair_count = int(fields[4]),
                junction_read_count = int(fields[5]),
                longest_anchor = int(fields[6]),
                methods = fields[7].split(';'),
                five_end_breakpoint = fields[8],
                three_end_breakpoint = fields[9],
                five_end_gene_id = fields[10],
                three_end_gene_id = fields[11],
                five_end_exon_id = fields[12],
                three_end_exon_id = fields[13],
                fusion_sequence = fields[14],
                predicted_effect = fields[15]
            )
            line = next(handle, None)


class FusionCatcherRecord():
    """ Defines a record of FusionCatcher's output. """
    def __init__(self, five_end_gene_symbol:str, three_end_gene_symbol:int,
            fusion_descriptions:List[str], multimapping_read_count:int,
            junction_pair_count:int, junction_read_count:int, longest_anchor: int,
            methods:List[str], five_end_breakpoint:str, three_end_breakpoint:str,
            five_end_gene_id:str, three_end_gene_id:str, five_end_exon_id:str,
            three_end_exon_id:str, fusion_sequence:str, predicted_effect:str
            ):
        """"""
        self.five_end_gene_symbol = five_end_gene_symbol
        self.three_end_gene_symbol = three_end_gene_symbol
        self.fusion_descriptions = fusion_descriptions
        self.multimapping_read_count = multimapping_read_count
        self.junction_pair_count = junction_pair_count
        self.junction_read_count = junction_read_count
        self.longest_anchor = longest_anchor
        self.methods = methods
        self.five_end_breakpoint = five_end_breakpoint
        self.three_end_breakpoint = three_end_breakpoint
        self.five_end_gene_id = five_end_gene_id
        self.three_end_gene_id = three_end_gene_id
        self.five_end_exon_id = five_end_exon_id
        self.three_end_exon_id = three_end_exon_id
        self.fusion_sequence = fusion_sequence
        self.predicted_effect = predicted_effect

    def convert_to_variant_records(self, anno:GeneAnnotations,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        """ COnvert a FusionCatcher's record to VariantRecord.

        Args:
            anno (Dict[str, Dict[str, gtf.TranscriptAnnotationModel]]): The
                genome annotation. The first level key is the gene ID, and the
                second level key is the transcript ID.
            genome (dna.DNASeqDict): The genome sequence records.
        Returns:
            List of VariantRecord
        """
        acceptor_transcripts = anno[self.three_end_gene_id]
        donor_transcripts = anno[self.five_end_gene_id]
        records = []

        perms = itertools.product(acceptor_transcripts.keys(), donor_transcripts.keys())
        for acceptor_id, donor_id in perms:
            # just in case ensembl ID does not match those in annotation
            if (acceptor_id is None or donor_id is None):
                continue
            # fusion catcher outputs results without 'chr' prefix, align with genome
            self.set_chrom_with_chr(genome.with_chr())

            acceptor_model = acceptor_transcripts[acceptor_id]
            acceptor_gene_symbol = acceptor_model.transcript.attributes['gene_name']
            # in case the ensembl ID does not match those in annotation
            # and not the same gene is referered to
            if (acceptor_gene_symbol != self.three_end_gene_symbol):
                continue
            # fusion catcher uses 1-based coordinates
            left_breakpoint = int(self.three_end_breakpoint.split(':')[1])
            left_strand = self.three_end_breakpoint.split(':')[2]
            acceptor_chrom = self.three_end_breakpoint.split(':')[0]
            acceptor_position = acceptor_model.get_transcript_index(left_breakpoint)
            seq = acceptor_model.get_transcript_sequence(genome[acceptor_chrom])
            acceptor_genome_position = f'{acceptor_chrom}:{left_breakpoint}:{left_breakpoint}'

            donor_model = donor_transcripts[donor_id]
            donor_gene_symbol = donor_model.transcript.attributes['gene_name']
            # in case the ensembl ID does not match those in annotation
            # and not the same gene is referered to
            if (donor_gene_symbol != self.five_end_gene_symbol):
                continue
            right_breakpoint = int(self.five_end_breakpoint.split(':')[1])
            right_strand = self.five_end_breakpoint.split(':')[2]

            donor_chrom = self.five_end_breakpoint.split(':')[0]
            donor_position = donor_model.get_transcript_index(right_breakpoint)
            donor_genome_position = f'{donor_chrom}:{right_breakpoint}:{right_breakpoint}'

            location = FeatureLocation(
                seqname=acceptor_id,
                start=acceptor_position,
                end=acceptor_position + 1
            )
            _id = f'FUSION_{acceptor_id}:{acceptor_position}-{donor_id}:{donor_position}'
            attrs = {
                'GENE_ID': self.three_end_gene_id,
                'GENE_SYMBOL': self.three_end_gene_symbol,
                'GENOMIC_POSITION': acceptor_genome_position,
                'DONOR_GENE_ID': self.five_end_gene_id,
                'DONOR_TRANSCRIPT_ID': donor_id,
                'DONOR_SYMBOL': self.five_end_gene_symbol,
                'DONOR_POS': donor_position,
                'DONOR_GENOMIC_POSITION': donor_genome_position
            }
            record = seqvar.VariantRecord(
                location=location,
                ref=seq[acceptor_position],
                alt='<FUSION>',
                _type='Fusion',
                _id=_id,
                attrs=attrs
            )
            records.append(record)
        return records

    def set_chrom_with_chr(self, with_chr=True) -> None:
        if (with_chr):
            if ('chr' not in self.five_end_breakpoint.split(':')[0]):
                self.five_end_breakpoint = 'chr' + self.five_end_breakpoint
                self.three_end_breakpoint = 'chr' + self.three_end_breakpoint
        if (not with_chr):
            if ('chr' in self.five_end_breakpoint.split(':')[0]):
                self.five_end_breakpoint[self.five_end_breakpoint.startswith('chr') and len('chr'):]
                self.three_end_breakpoint[self.three_end_breakpoint.startswith('chr') and len('chr'):]
