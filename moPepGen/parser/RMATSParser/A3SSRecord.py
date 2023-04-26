""" Module for rMATS parser """
from __future__ import annotations
from typing import TYPE_CHECKING, List, Tuple
from moPepGen import gtf, dna, seqvar
from .RMATSRecord import RMATSRecord


if TYPE_CHECKING:
    from moPepGen.seqvar import SpliceJunction


class A3SSRecord(RMATSRecord):
    """ Alternative 3' Splicing Site """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            long_exon_start:int, long_exon_end:int, short_exon_start:int,
            short_exon_end:int, flanking_exon_start, flanking_exon_end:int,
            ijc_sample_1:int, sjc_sample_1:int, ijc_sample_2:int,
            sjc_sample_2:int, inc_form_len:int, skip_form_len:int,
            pvalue:float, fdr:float):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.long_exon_start = long_exon_start
        self.long_exon_end = long_exon_end
        self.short_exon_start = short_exon_start
        self.short_exon_end = short_exon_end
        self.flanking_exon_start = flanking_exon_start
        self.flanking_exon_end = flanking_exon_end
        self.ijc_sample_1 = ijc_sample_1
        self.sjc_sample_1 = sjc_sample_1
        self.ijc_sample_2 = ijc_sample_2
        self.sjc_sample_2 = sjc_sample_2
        self.inc_form_len = inc_form_len
        self.skip_form_len = skip_form_len
        self.pvalue = pvalue
        self.fdr = fdr

    @classmethod
    def readline(cls, line:str) -> A3SSRecord:
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
            flanking_exon_end=int(fields[10]),
            ijc_sample_1=int(fields[12]),
            sjc_sample_1=int(fields[13]),
            ijc_sample_2=None if fields[14] == '' else int(fields[14]),
            sjc_sample_2=None if fields[15] == '' else int(fields[15]),
            inc_form_len=int(fields[16]),
            skip_form_len=int(fields[17]),
            pvalue=None if fields[18] == 'NA' else float(fields[18]),
            fdr=None if fields[19] == 'NA' else float(fields[19])
        )

    def create_variant_id(self, anno:gtf.GenomicAnnotation) -> str:
        """ Create variant ID """
        if anno.genes[self.gene_id].strand == 1:
            uee = anno.coordinate_genomic_to_gene(self.flanking_exon_end - 1, self.gene_id) + 1
            les = anno.coordinate_genomic_to_gene(self.long_exon_start, self.gene_id)
            ses = anno.coordinate_genomic_to_gene(self.short_exon_start, self.gene_id)
        else:
            uee = anno.coordinate_genomic_to_gene(self.flanking_exon_start, self.gene_id) + 1
            les = anno.coordinate_genomic_to_gene(self.long_exon_end - 1, self.gene_id)
            ses = anno.coordinate_genomic_to_gene(self.short_exon_end - 1, self.gene_id)

        return f"A5SS_{uee}-{les}-{ses}"

    def create_splice_junctions(self, strand:int) -> Tuple[SpliceJunction, SpliceJunction]:
        """ Create splice junctions """
        if strand == 1:
            long_junction = seqvar.SpliceJunction(
                upstream_start=self.flanking_exon_start,
                upstream_end=self.flanking_exon_end,
                downstream_start=self.long_exon_start,
                downstream_end=self.long_exon_end,
                gene_id=self.gene_id,
                chrom=self.chrom
            )
            short_junction = seqvar.SpliceJunction(
                upstream_start=self.flanking_exon_start,
                upstream_end=self.flanking_exon_end,
                downstream_start=self.short_exon_start,
                downstream_end=self.short_exon_end,
                gene_id=self.gene_id,
                chrom=self.chrom
            )
        else:
            long_junction = seqvar.SpliceJunction(
                upstream_start=self.long_exon_start,
                upstream_end=self.long_exon_end,
                downstream_start=self.flanking_exon_start,
                downstream_end=self.flanking_exon_end,
                gene_id=self.gene_id,
                chrom=self.chrom
            )
            short_junction = seqvar.SpliceJunction(
                upstream_start=self.short_exon_start,
                upstream_end=self.short_exon_end,
                downstream_start=self.flanking_exon_start,
                downstream_end=self.flanking_exon_end,
                gene_id=self.gene_id,
                chrom=self.chrom
            )
        return long_junction, short_junction

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        """ Convert A3SS record to VariantRecord """
        variants = []

        gene_model = anno.genes[self.gene_id]
        strand = gene_model.strand

        long_junction, short_junction = self.create_splice_junctions(strand)

        if not long_junction.is_novel(anno) and not short_junction.is_novel(anno):
            return variants

        tx_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])
        var_id = self.create_variant_id(anno)

        if gene_model.strand == 1:
            upstream_novel = False
            downstream_novel = True
        else:
            upstream_novel = True
            downstream_novel = False

        # For A3SS, long is inclusion
        for tx_id in tx_ids:
            tx_model = anno.transcripts[tx_id]

            if self.ijc_sample_1 >= min_ijc:
                aln = long_junction.align_to_transcript(
                    tx_model, upstream_novel, downstream_novel
                )
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)

            if self.sjc_sample_1 >= min_sjc:
                aln = short_junction.align_to_transcript(
                    tx_model, upstream_novel, downstream_novel
                )
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)
        return variants
