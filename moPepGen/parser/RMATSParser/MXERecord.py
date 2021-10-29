""" Module fro Mutually Exclusive Exons """
from __future__ import annotations
from typing import List
from moPepGen import seqvar, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord

class MXERecord(RMATSRecord):
    """ Mutually Exclusive Exons """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            first_exon_start:int, first_exon_end:int, second_exon_start:int,
            second_exon_end:int, upstream_exon_start:int,
            upstream_exon_end:int, downstream_exon_start:int,
            downstream_exon_end:int):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.first_exon_start = first_exon_start
        self.first_exon_end = first_exon_end
        self.second_exon_start = second_exon_start
        self.second_exon_end = second_exon_end
        self.upstream_exon_start = upstream_exon_start
        self.upstream_exon_end = upstream_exon_end
        self.downstream_exon_start = downstream_exon_start
        self.downstream_exon_end = downstream_exon_end

    @classmethod
    def readline(cls, line:str) -> MXERecord:
        """ Create a A5SSRecord from a line of the rMASTS output """
        fields = line.rstrip().split('\t')
        return cls(
            gene_id=fields[1].strip('"'),
            gene_symbol=fields[2].strip('"'),
            chrom=fields[3],
            first_exon_start=int(fields[5]),
            first_exon_end=int(fields[6]),
            second_exon_start=int(fields[7]),
            second_exon_end=int(fields[8]),
            upstream_exon_start=int(fields[9]),
            upstream_exon_end=int(fields[10]),
            downstream_exon_start=int(fields[11]),
            downstream_exon_end=int(fields[12])
        )

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        variants = []
        gene_model = anno.genes[self.gene_id]
        tx_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        strand = gene_model.location.strand
        gene_seq = gene_model.get_gene_sequence(genome[chrom])

        have_both:List[str] = []
        have_first:List[str] = []
        have_second:List[str] = []

        for tx_id in tx_ids:
            model = anno.transcripts[tx_id]
            it = iter(model.exon)
            exon = next(it, None)
            while exon:
                if int(exon.location.end) != self.upstream_exon_end:
                    exon = next(it, None)
                    continue
                exon = next(it, None)
                if not exon:
                    break
                if int(exon.location.start) == self.first_exon_start and \
                        int(exon.location.end) == self.first_exon_end:
                    exon = next(it, None)
                    if not exon:
                        break
                    if int(exon.location.start) == self.downstream_exon_start:
                        if strand == 1:
                            have_first.append(tx_id)
                        else:
                            have_second.append(tx_id)
                    elif int(exon.location.start) == self.second_exon_start and \
                            int(exon.location.end) == self.second_exon_end:
                        exon = next(it, None)
                        if not exon:
                            break
                        if int(exon.location.start) == self.downstream_exon_start:
                            have_both.append(tx_id)
                elif int(exon.location.start) == self.second_exon_start and \
                        int(exon.location.end) == self.second_exon_end:
                    exon = next(it, None)
                    if not exon:
                        break
                    if int(exon.location.start) == self.downstream_exon_start:
                        if strand == 1:
                            have_second.append(tx_id)
                        else:
                            have_first.append(tx_id)

        if (have_first and have_second) or \
                (not have_first and not have_second and not have_both):
            return variants

        if strand == 1:
            first_start = anno.coordinate_genomic_to_gene(
                self.first_exon_start, self.gene_id)
            first_end = anno.coordinate_genomic_to_gene(
                self.first_exon_end - 1, self.gene_id) + 1
            second_start = anno.coordinate_genomic_to_gene(
                self.second_exon_start, self.gene_id)
            second_end = anno.coordinate_genomic_to_gene(
                self.second_exon_end - 1, self.gene_id) + 1
            first_genomic_position = f'{chrom}:{self.first_exon_start + 1}'+\
                f'-{self.first_exon_end}'
            second_genomic_position = f'{chrom}:{self.second_exon_start + 1}'+\
                f'-{self.second_exon_end}'
        else:
            first_start = anno.coordinate_genomic_to_gene(
                self.second_exon_end - 1, self.gene_id)
            first_end = anno.coordinate_genomic_to_gene(
                self.second_exon_start, self.gene_id) + 1
            second_start = anno.coordinate_genomic_to_gene(
                self.first_exon_end - 1, self.gene_id)
            second_end = anno.coordinate_genomic_to_gene(
                self.first_exon_start, self.gene_id) + 1
            second_genomic_position = f'{chrom}:{self.first_exon_start + 1}'+\
                f'-{self.first_exon_end}'
            first_genomic_position = f'{chrom}:{self.second_exon_end + 1}'+\
                f'-{self.second_exon_end}'

        _id = f'MXE_{first_start + 1}-{first_end}:' +\
            f'{second_start}-{second_end}'

        if not have_second:
            for tx_id in have_first:
                location = FeatureLocation(seqname=self.gene_id, start=first_start,
                    end=first_end)
                ref = str(gene_seq.seq[first_start])
                alt = '<SUB>'
                attrs = {
                    'GENE_ID': tx_id,
                    'START': first_start,
                    'END': first_end,
                    'DONOR_START': second_start,
                    'DONOR_END': second_end,
                    'DONOR_GENE_ID': self.gene_id,
                    'COORDINATE': 'gene',
                    'GENE_SYMBOL': model.transcript.gene_name,
                    'GENOMIC_POSITION': f'{chrom}-{first_start + 1}:{first_end}-'
                    f'{second_start + 1}:{second_end}'
                }
                _type = 'Substitution'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

            for tx_id in have_both:
                location = FeatureLocation(seqname=self.gene_id, start=first_start,
                    end=first_end)
                ref = str(gene_seq.seq[first_start])
                alt = '<DEL>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'START': first_start,
                    'END': first_end,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': first_genomic_position
                }
                _type = 'Deletion'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not have_first:
            for tx_id in have_second:
                location = FeatureLocation(seqname=self.gene_id, start=second_start,
                    end=second_end)
                ref = str(gene_seq.seq[second_start])
                alt = '<SUB>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'START': second_start,
                    'END': second_end,
                    'DONOR_START': first_start,
                    'DONOR_END': first_end,
                    'DONOR_GENE_ID': self.gene_id,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': f'{chrom}-{second_start + 1}:{second_end}-'
                    f'{first_start + 1}:{first_end}'
                }
                _type = 'Substitution'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

            for tx_id in have_both:
                location = FeatureLocation(seqname=self.gene_id, start=second_start,
                    end=second_end)
                ref = str(gene_seq.seq[second_start])
                alt = '<DEL>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'START': second_start,
                    'END': second_end,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': second_genomic_position
                }
                _type = 'Deletion'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        return variants
