""" Module for FusionCatcher parser """
from __future__ import annotations
import re
from typing import Iterable, List, Tuple
import itertools
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, gtf, dna, err


def parse(path:str) -> Iterable[FusionCatcherRecord]:
    """ Parse the FusionCatcher's output and returns an iterator.

    Args:
        path (str): Path to the FusionCatcher's output file.
    """
    with open(path, 'r') as handle:
        # first line is header
        line = next(handle, None)
        line = next(handle, None)
        while line:
            if line.startswith('#'):
                line = next(handle, None)
                continue
            fields = line.rstrip().split('\t')
            yield FusionCatcherRecord(
                five_end_gene_symbol=fields[0],
                three_end_gene_symbol=fields[1],
                fusion_descriptions=fields[2].split(','),
                counts_of_common_mapping_reads=int(fields[3]),
                spanning_pairs=int(fields[4]),
                spanning_unique_reads=int(fields[5]),
                longest_anchor_found=int(fields[6]),
                fusion_finding_method=fields[7].split(';'),
                five_end_breakpoint=fields[8],
                three_end_breakpoint=fields[9],
                five_end_gene_id=fields[10],
                three_end_gene_id=fields[11],
                five_end_exon_id=fields[12],
                three_end_exon_id=fields[13],
                fusion_sequence=tuple(fields[14].split('*')),
                predicted_effect=fields[15]
            )
            line = next(handle, None)


class FusionCatcherRecord():
    """ Defines a record of FusionCatcher's output. """
    def __init__(self, five_end_gene_symbol:str, three_end_gene_symbol:str,
            fusion_descriptions:List[str], counts_of_common_mapping_reads:int,
            spanning_pairs:int, spanning_unique_reads:int,
            longest_anchor_found:int, fusion_finding_method:List[str],
            five_end_breakpoint:str, three_end_breakpoint:str,
            five_end_gene_id:str, three_end_gene_id:str, five_end_exon_id:str,
            three_end_exon_id:str, fusion_sequence:Tuple[str, str],
            predicted_effect:str
            ):
        """"""
        self.five_end_gene_symbol = five_end_gene_symbol
        self.three_end_gene_symbol = three_end_gene_symbol
        self.fusion_descriptions = fusion_descriptions
        self.counts_of_common_mapping_reads = counts_of_common_mapping_reads
        self.spanning_pairs = spanning_pairs
        self.spanning_unique_reads = spanning_unique_reads
        self.longest_anchor_found = longest_anchor_found
        self.fusion_finding_method = fusion_finding_method
        self.five_end_breakpoint = five_end_breakpoint
        self.three_end_breakpoint = three_end_breakpoint
        self.five_end_gene_id = five_end_gene_id
        self.three_end_gene_id = three_end_gene_id
        self.five_end_exon_id = five_end_exon_id
        self.three_end_exon_id = three_end_exon_id
        self.fusion_sequence = fusion_sequence
        self.predicted_effect = predicted_effect

    @property
    def left_breakpoint_position(self) -> int:
        """ Get left breakpoint position """
        return int(self.five_end_breakpoint.split(':')[1])

    @property
    def right_breakpoint_position(self) -> int:
        """ Get right breakpoint position """
        return int(self.three_end_breakpoint.split(':')[1])

    def get_donor_transcripts(self, anno:gtf.GenomicAnnotation, gene_id:str
            ) -> List[gtf.TranscriptAnnotationModel]:
        """ Get all possible donor transcripts """
        pos = self.left_breakpoint_position - 1
        return anno.get_transcripts_with_position(gene_id, pos)

    def get_accepter_transcripts(self, anno:gtf.GenomicAnnotation, gene_id:str
            ) -> List[gtf.TranscriptAnnotationModel]:
        """ Get all possible accepter transcripts """
        pos = self.right_breakpoint_position - 1
        return anno.get_transcripts_with_position(gene_id, pos)

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
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
        pattern = re.compile(r'\.[0-9]+$')
        if pattern.search(self.five_end_gene_id):
            donor_gene_id = self.five_end_gene_id
            accepter_gene_id = self.three_end_gene_id
            try:
                donor_gene_model = anno.genes[self.five_end_gene_id]
            except KeyError as error:
                raise err.GeneNotFoundError(self.five_end_gene_id) from error
            try:
                accepter_gene_model = anno.genes[self.three_end_gene_id]
            except KeyError as error:
                raise err.GeneNotFoundError(self.three_end_gene_id) from error
        else:
            donor_gene_model = anno.get_gene_model_from_unversioned_id(
                self.five_end_gene_id)
            accepter_gene_model = anno.get_gene_model_from_unversioned_id(
                self.three_end_gene_id)
            donor_gene_id = donor_gene_model.gene_id
            accepter_gene_id = accepter_gene_model.gene_id

        # in case the ensembl ID does not match those in annotation
        # and not the same gene is referered to
        donor_gene_symbol = donor_gene_model.gene_name
        accepter_gene_symbol = accepter_gene_model.gene_name

        donor_chrom = donor_gene_model.chrom
        accepter_chrom = accepter_gene_model.chrom

        # fusion catcher uses 1-based coordinates
        # left breakpoint is the first nucleotide in the fusion transcript
        # after the breakpoint
        left_breakpoint_genomic = self.left_breakpoint_position
        right_breakpoint_genomic = self.right_breakpoint_position
        donor_genome_position = \
            f'{donor_chrom}:{left_breakpoint_genomic}:{left_breakpoint_genomic}'
        accepter_genome_position = \
            f'{accepter_chrom}:{right_breakpoint_genomic}:{right_breakpoint_genomic}'
        left_breakpoint_genetic = anno.coordinate_genomic_to_gene(
            index=left_breakpoint_genomic - 1, gene=donor_gene_id
        ) + 1
        right_breakpoint_genetic = anno.coordinate_genomic_to_gene(
            index=right_breakpoint_genomic - 1, gene=accepter_gene_id
        )

        donor_transcripts = self.get_donor_transcripts(anno, donor_gene_id)
        accepter_transcripts = self.get_accepter_transcripts(anno, accepter_gene_id)

        records = []

        if donor_gene_model.strand == 1:
            ref_seq = genome[donor_chrom].seq[left_breakpoint_genomic]
        else:
            ref_seq = genome[donor_chrom]\
                .seq[left_breakpoint_genomic:left_breakpoint_genomic + 1]\
                .reverse_complement()
            ref_seq = str(ref_seq)

        perms = itertools.product(donor_transcripts, accepter_transcripts)
        for donor_tx, accepter_tx in perms:
            donor_tx_id = donor_tx.transcript.transcript_id
            accepter_tx_id = accepter_tx.transcript.transcript_id

            fusion_id = f'FUSION-{donor_tx_id}:{left_breakpoint_genetic}'+\
                f'-{accepter_tx_id}:{right_breakpoint_genetic}'

            location = FeatureLocation(
                seqname=donor_gene_id,
                start=left_breakpoint_genetic,
                end=left_breakpoint_genetic + 1
            )
            attrs = {
                'TRANSCRIPT_ID': donor_tx_id,
                'GENE_SYMBOL': donor_gene_symbol,
                'GENOMIC_POSITION': donor_genome_position,
                'ACCEPTER_GENE_ID': accepter_gene_id,
                'ACCEPTER_TRANSCRIPT_ID': accepter_tx_id,
                'ACCEPTER_SYMBOL': accepter_gene_symbol,
                'ACCEPTER_POSITION': right_breakpoint_genetic,
                'ACCEPTER_GENOMIC_POSITION': accepter_genome_position
            }

            record = seqvar.VariantRecord(
                location=location,
                ref=ref_seq,
                alt='<FUSION>',
                _type='Fusion',
                _id=fusion_id,
                attrs=attrs
            )
            records.append(record)
        return records

    def set_chrom_with_chr(self, with_chr=True) -> None:
        """ Depending on if chr is required in chrom, fix breakpoints"""
        if with_chr:
            if 'chr' not in self.five_end_breakpoint.split(':')[0]:
                self.five_end_breakpoint = 'chr' + self.five_end_breakpoint
                self.three_end_breakpoint = 'chr' + self.three_end_breakpoint
        if not with_chr:
            if 'chr' in self.five_end_breakpoint.split(':')[0]:
                feb = self.five_end_breakpoint
                feb = feb[feb.startswith('chr') and len('chr'):]
                teb = self.three_end_breakpoint
                teb = teb[teb.startswith('chr') and len('chr'):]
                self.five_end_breakpoint = feb
                self.three_end_breakpoint = teb
