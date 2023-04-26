""" Module for rMATS parser """
from __future__ import annotations
from typing import TYPE_CHECKING, List, Tuple
from moPepGen import gtf, dna, seqvar
from .RMATSRecord import RMATSRecord


if TYPE_CHECKING:
    from moPepGen.seqvar import SpliceJunction


class SERecord(RMATSRecord):
    """ Skipped Exon """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str, exon_start:int,
            exon_end:int, upstream_exon_start:int, upstream_exon_end:int,
            downstream_exon_start:int, downstream_exon_end:int,
            ijc_sample_1:int, sjc_sample_1:int, ijc_sample_2:int,
            sjc_sample_2:int, inc_form_len:int, skip_form_len:int,
            pvalue:float, fdr:float):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.upstream_exon_start = upstream_exon_start
        self.upstream_exon_end = upstream_exon_end
        self.downstream_exon_start = downstream_exon_start
        self.downstream_exon_end = downstream_exon_end
        self.ijc_sample_1 = ijc_sample_1
        self.sjc_sample_1 = sjc_sample_1
        self.ijc_sample_2 = ijc_sample_2
        self.sjc_sample_2 = sjc_sample_2
        self.inc_form_len = inc_form_len
        self.skip_form_len = skip_form_len
        self.pvalue = pvalue
        self.fdr = fdr

    @classmethod
    def readline(cls, line:str) -> SERecord:
        """ Create a SERecord from a line of the rMATS output """
        fields = line.rstrip().split('\t')
        return cls(
            gene_id=fields[1].strip('"'),
            gene_symbol=fields[2].strip('"'),
            chrom=fields[3],
            exon_start=int(fields[5]),
            exon_end=int(fields[6]),
            upstream_exon_start=int(fields[7]),
            upstream_exon_end=int(fields[8]),
            downstream_exon_start=int(fields[9]),
            downstream_exon_end=int(fields[10]),
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
        uee = anno.coordinate_genomic_to_gene(self.upstream_exon_end - 1, self.gene_id)
        es = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        ee = anno.coordinate_genomic_to_gene(self.exon_end - 1, self.gene_id)
        des = anno.coordinate_genomic_to_gene(self.downstream_exon_start, self.gene_id)

        gene_model = anno.genes[self.gene_id]
        if gene_model.strand == -1:
            es, ee = ee, es
            uee, des = des, uee
        ee += 1
        des += 1

        return f"SE_{uee}-{es}-{ee}-{des}"

    def create_splice_junctions(self
            ) -> Tuple[SpliceJunction, SpliceJunction, SpliceJunction]:
        """ Create splice junctions """
        skip_junction = seqvar.SpliceJunction(
            upstream_start=self.upstream_exon_start,
            upstream_end=self.upstream_exon_end,
            downstream_start=self.downstream_exon_start,
            downstream_end=self.downstream_exon_end,
            gene_id=self.gene_id,
            chrom=self.chrom
        )
        upstream_junction = seqvar.SpliceJunction(
            upstream_start=self.upstream_exon_start,
            upstream_end=self.upstream_exon_end,
            downstream_start=self.exon_start,
            downstream_end=self.exon_end,
            gene_id=self.gene_id,
            chrom=self.chrom
        )
        downstream_junction = seqvar.SpliceJunction(
            upstream_start=self.exon_start,
            upstream_end=self.exon_end,
            downstream_start=self.downstream_exon_start,
            downstream_end=self.downstream_exon_end,
            gene_id=self.gene_id,
            chrom=self.chrom
        )
        return skip_junction, upstream_junction, downstream_junction

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        """ Convert to variant records """
        variants = []

        gene_model = anno.genes[self.gene_id]

        skip_junction, upstream_junction, downstream_junction \
            = self.create_splice_junctions()

        if not skip_junction.is_novel(anno) \
                and not upstream_junction.is_novel(anno)\
                and not downstream_junction.is_novel(anno):
            return variants

        tx_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])
        var_id = self.create_variant_id(anno)

        for tx_id in tx_ids:
            tx_model = anno.transcripts[tx_id]

            if self.sjc_sample_1 >= min_sjc:
                aln = skip_junction.align_to_transcript(tx_model, False, False)
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)

            if self.ijc_sample_1 >= min_ijc:
                aln = upstream_junction.align_to_transcript( tx_model, False, True)
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)
                aln = downstream_junction.align_to_transcript(tx_model, True, False)
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)

        return variants
