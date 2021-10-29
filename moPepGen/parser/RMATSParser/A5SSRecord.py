""" Module for Altervative 5' Splicing Site """
from __future__ import annotations
from typing import List
from moPepGen import seqvar, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord


class A5SSRecord(RMATSRecord):
    """ Alternative 5' Splicing Site """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            long_exon_start:int, long_exon_end:int, short_exon_start:int,
            short_exon_end:int, flanking_exon_start, flanking_exon_end:int):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.long_exon_start = long_exon_start
        self.long_exon_end = long_exon_end
        self.short_exon_start = short_exon_start
        self.short_exon_end = short_exon_end
        self.flanking_exon_start = flanking_exon_start
        self.flanking_exon_end = flanking_exon_end

    @classmethod
    def readline(cls, line:str) -> A5SSRecord:
        """ Create a A5SSRecord from a line of the rMASTS output """
        fields = line.rstrip().split('\t')
        return cls(
            gene_id=fields[1].strip('"'),
            gene_symbol=fields[2].strip('"'),
            chrom=fields[3],
            long_exon_start=int(fields[5]),
            long_exon_end=int(fields[6]),
            short_exon_start=int(fields[7]),
            short_exon_end=int(fields[8]),
            flanking_exon_start=int(fields[9]),
            flanking_exon_end=int(fields[10])
        )

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        variants = []
        gene_model = anno.genes[self.gene_id]
        tx_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])

        short:List[str] = []
        long:List[str] = []

        for tx_id in tx_ids:
            model = anno.transcripts[tx_id]
            if gene_model.location.strand == 1:
                it = iter(model.exon)
            else:
                it = reversed(model.exon)
            exon = next(it, None)

            while exon:
                if int(exon.location.start) == self.long_exon_start and \
                    int(exon.location.end) == self.long_exon_end:
                    exon = next(it, None)
                    if exon and exon.location.start == self.flanking_exon_start:
                        long.append(tx_id)
                    break
                if int(exon.location.start) == self.short_exon_start and \
                        int(exon.location.end) == self.short_exon_end:
                    exon = next(it, None)
                    if exon and exon.location.start == self.flanking_exon_start:
                        short.append(tx_id)
                    break
                exon = next(it, None)

        if (long and short) or (not long and not short):
            return variants

        if anno.genes[self.gene_id].location.strand == 1:
            start_gene = anno.coordinate_genomic_to_gene(self.short_exon_end - 1,
                self.gene_id)
            end_gene = anno.coordinate_genomic_to_gene(self.long_exon_end - 1,
                self.gene_id)
            genomic_position = f'{chrom}:{self.short_exon_end+1}-{self.long_exon_end}'
        else:
            start_gene = anno.coordinate_genomic_to_gene(self.short_exon_start,
                self.gene_id)
            end_gene = anno.coordinate_genomic_to_gene(self.long_exon_start,
                self.gene_id)
            genomic_position = f'{chrom}:{self.long_exon_start+1}-{self.short_exon_end}'

        if not short:
            for tx_id in long:
                location = FeatureLocation(seqname=self.gene_id, start=start_gene,
                    end=end_gene)
                ref = str(gene_seq.seq[start_gene])
                alt = '<DEL>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'START': start_gene,
                    'END': end_gene,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Deletion'
                _id = f'A5SS_{start_gene}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not long:
            insert_position = start_gene
            for tx_id in short:
                location = FeatureLocation(
                    seqname=self.gene_id,
                    start=insert_position,
                    end=insert_position+1
                )
                ref = str(gene_seq.seq[insert_position])
                alt = '<INS>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'DONOR_START': start_gene,
                    'DONOR_END': end_gene,
                    'DONOR_GENE_ID': self.gene_id,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Insertion'
                _id = f'SE_{start_gene}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        return variants
