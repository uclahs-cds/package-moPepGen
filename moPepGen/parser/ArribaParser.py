""" Module for Arriba parser """
from __future__ import annotations
import itertools
from pathlib import Path
from typing import Iterable, List, IO
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, gtf, dna, err


def parse(handle:IO) -> Iterable[ArribaRecord]:
    """ parse Arriba's output """
    for line in handle:
        line:str
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        fields[9] = int(fields[9])
        fields[10] = int(fields[10])
        fields[11] = int(fields[11])
        fields[12] = int(fields[12])
        fields[13] = int(fields[13])
        fields[26] = fields[26].split(',')
        fields[29] = fields[29].split(',')

        yield ArribaRecord(*fields)

class ArribaRecord():
    """ Arriba Record """
    def __init__(self, gene1:str, gene2:str, strand1:str, strand2:str,
            breakpoint1:str, breakpoint2:str, site1:str, site2:str, type:str,
            split_reads1:str, split_reads2:str, discordant_mates:int,
            coverage1:int, coverage2:int, confidence:str, reading_frame:str,
            tags:str, retained_protein_domains:str,
            closest_genomic_breakpoint1:str, closest_genomic_breakpoint2:str,
            gene_id1:str, gene_id2:str, transcript_id1:str, transcript_id2:str,
            direction1:str, direction2:str, filters:List[str],
            fusion_transcript:str, peptide_sequence:str,
            read_identifiers:List[str]):
        """ Constructor """
        self.gene1 = gene1
        self.gene2 = gene2
        self.strand1 = strand1
        self.strand2 = strand2
        self.breakpoint1 = breakpoint1
        self.breakpoint2 = breakpoint2
        self.site1 = site1
        self.site2 = site2
        self.type = type
        self.split_reads1 = split_reads1
        self.split_reads2 = split_reads2
        self.discordant_mates = discordant_mates
        self.coverage1 = coverage1
        self.coverage2 = coverage2
        self.confidence = confidence
        self.reading_frame = reading_frame
        self.tags = tags
        self.retained_protein_domains = retained_protein_domains
        self.closest_genomic_breakpoint1 = closest_genomic_breakpoint1
        self.closest_genomic_breakpoint2 = closest_genomic_breakpoint2
        self.gene_id1 = gene_id1
        self.gene_id2 = gene_id2
        self.transcript_id1 = transcript_id1
        self.transcript_id2 = transcript_id2
        self.direction1 = direction1
        self.direction2 = direction2
        self.filters = filters
        self.fusion_transcript = fusion_transcript
        self.peptide_sequence = peptide_sequence
        self.read_identifiers = read_identifiers

    @property
    def gene_strand1(self) -> int:
        """ gene strand1 """
        return self.infer_strand(self.strand1.split('/')[0])

    @property
    def gene_strand2(self) -> int:
        """ gene strand2 """
        return self.infer_strand(self.strand2.split('/')[0])

    @property
    def transcript_strand1(self) -> int:
        """ transcript strand1 """
        return self.infer_strand(self.strand1.split('/')[1])

    @property
    def transcript_strand2(self) -> int:
        """ transcript strand2 """
        return self.infer_strand(self.strand2.split('/')[1])

    @property
    def breakpoint1_position(self) -> int:
        """ Get breakpoint1 position """
        return int(self.breakpoint1.split(':')[1])

    @property
    def breakpoint2_position(self) -> int:
        """ Get breakpoint2 position """
        return int(self.breakpoint2.split(':')[1])

    @staticmethod
    def infer_strand(strand:str):
        """ Infer strand """
        if strand == '+':
            return 1
        if strand == '-':
            return -1
        if strand == '.':
            return 0
        raise ValueError(f"Can not infer strand from {strand}")

    def transcript_on_antisense_strand(self, anno:gtf.GenomicAnnotation) -> bool:
        """ """
        gene_model1 = anno.genes[self.gene_id1]
        gene_model2 = anno.genes[self.gene_id2]
        return self.transcript_strand1 != gene_model1.strand or\
            self.transcript_strand2 != gene_model2.strand

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        """ """
        try:
            donor_gene_model = anno.genes[self.gene1]
        except KeyError as error:
            raise err.GeneNotFoundError(self.gene1) from error
        try:
            accepter_gene_model = anno.genes[self.gene2]
        except KeyError as error:
            raise err.GeneNotFoundError(self.gene2) from error

        donor_gene_symbol = donor_gene_model.gene_name
        donor_chrom = donor_gene_model.chrom
        left_breakpoint = self.breakpoint1_position
        donor_genome_position = f'{donor_chrom}:{left_breakpoint}:{left_breakpoint}'
        donor_position = anno.coordinate_genomic_to_gene(left_breakpoint - 1, self.gene1) + 1
        donor_transcripts = [x for x in donor_gene_model.transcripts
             if anno.transcripts[x].is_exonic(left_breakpoint - 1)]

        accepter_gene_symbol = accepter_gene_model.gene_name
        accepter_chrom = accepter_gene_model.chrom
        right_breakpoint = self.breakpoint2_position
        accepter_genome_position = f'{accepter_chrom}:{right_breakpoint}:{right_breakpoint}'
        accepter_position = anno.coordinate_genomic_to_gene(right_breakpoint - 1, self.gene1)
        accepter_transcripts = [x for x in accepter_gene_model.transcripts
            if anno.transcripts[x].is_exonic(right_breakpoint - 1)]

        records = []

        fusion_id = f'FUSION-{self.left_gene}:{donor_position}'\
            f'-{self.right_gene}:{accepter_position}'

        if donor_gene_model.strand == 1:
            ref_seq = genome[donor_chrom].seq[donor_genome_position + 1]
        else:
            ref_seq = genome[donor_chrom]\
                .seq[donor_genome_position - 1:donor_genome_position]\
                .reverse_complement()
            ref_seq = str(ref_seq)

        perms = itertools.product(donor_transcripts, accepter_transcripts)
        for donor_tx, accepter_tx in perms:
            location = FeatureLocation(
                seqname=self.gene_id1,
                start=donor_position,
                end=donor_position + 1
            )
            attrs = {
                'TRANSCRIPT_ID': donor_tx,
                'GENE_SYMBOL': donor_gene_symbol,
                'GENOMIC_POSITION': donor_genome_position,
                'ACCEPTER_GENE_ID': self.right_gene,
                'ACCEPTER_TRANSCRIPT_ID': accepter_tx,
                'ACCEPTER_SYMBOL': accepter_gene_symbol,
                'ACCEPTER_POSITION': accepter_position,
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
