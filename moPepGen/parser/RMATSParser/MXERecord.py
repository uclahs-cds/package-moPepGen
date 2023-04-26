""" Module fro Mutually Exclusive Exons """
from __future__ import annotations
from typing import TYPE_CHECKING, List, Tuple
from moPepGen import gtf, dna, seqvar
from .RMATSRecord import RMATSRecord


if TYPE_CHECKING:
    from moPepGen.seqvar import SpliceJunction


class MXERecord(RMATSRecord):
    """ Mutually Exclusive Exons """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            first_exon_start:int, first_exon_end:int, second_exon_start:int,
            second_exon_end:int, upstream_exon_start:int,
            upstream_exon_end:int, downstream_exon_start:int,
            downstream_exon_end:int, ijc_sample_1:int, sjc_sample_1:int,
            ijc_sample_2:int, sjc_sample_2:int, inc_form_len:int,
            skip_form_len:int, pvalue:float, fdr:float):
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
        self.ijc_sample_1 = ijc_sample_1
        self.sjc_sample_1 = sjc_sample_1
        self.ijc_sample_2 = ijc_sample_2
        self.sjc_sample_2 = sjc_sample_2
        self.inc_form_len = inc_form_len
        self.skip_form_len = skip_form_len
        self.pvalue = pvalue
        self.fdr = fdr

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
            downstream_exon_end=int(fields[12]),
            ijc_sample_1=int(fields[14]),
            sjc_sample_1=int(fields[15]),
            ijc_sample_2=None if fields[16] == '' else int(fields[16]),
            sjc_sample_2=None if fields[17] == '' else int(fields[17]),
            inc_form_len=int(fields[18]),
            skip_form_len=int(fields[19]),
            pvalue=None if fields[20] == 'NA' else float(fields[20]),
            fdr=None if fields[21] == 'NA' else float(fields[21])
        )

    def create_variant_id(self, anno:gtf.GenomicAnnotation) -> str:
        """ Create variant ID """
        uee = anno.coordinate_genomic_to_gene(self.upstream_exon_end - 1, self.gene_id)
        fes = anno.coordinate_genomic_to_gene(self.first_exon_start, self.gene_id)
        fee = anno.coordinate_genomic_to_gene(self.first_exon_end - 1, self.gene_id)
        ses = anno.coordinate_genomic_to_gene(self.second_exon_start, self.gene_id)
        see = anno.coordinate_genomic_to_gene(self.second_exon_end - 1, self.gene_id)
        des = anno.coordinate_genomic_to_gene(self.downstream_exon_start, self.gene_id)

        gene_model = anno.genes[self.gene_id]
        if gene_model.strand == -1:
            uee, des = des, uee
            fes, fee = fee, fes
            ses, see = see, ses

        uee += 1
        fee += 1
        see += 1

        return f"MXE_{uee}-{fes}-{fee}-{ses}-{see}-{des}"

    def create_splice_junctions(self
            ) -> Tuple[SpliceJunction, SpliceJunction, SpliceJunction, SpliceJunction]:
        """ Create splice junctions """
        first_downstream_junction = seqvar.SpliceJunction(
            upstream_start=self.first_exon_start,
            upstream_end=self.first_exon_end,
            downstream_start=self.downstream_exon_start,
            downstream_end=self.downstream_exon_end,
            gene_id=self.gene_id,
            chrom=self.chrom
        )
        second_upstream_junction = seqvar.SpliceJunction(
            upstream_start=self.upstream_exon_start,
            upstream_end=self.upstream_exon_end,
            downstream_start=self.second_exon_start,
            downstream_end=self.second_exon_end,
            gene_id=self.gene_id,
            chrom=self.chrom
        )
        return first_downstream_junction, second_upstream_junction

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        variants = []

        gene_model = anno.genes[self.gene_id]

        first_downstream_junction, second_upstream_junction \
            = self.create_splice_junctions()

        if not first_downstream_junction.is_novel(anno)\
                and not second_upstream_junction.is_novel(anno):
            return variants

        tx_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])
        var_id = self.create_variant_id(anno)

        for tx_id in tx_ids:
            tx_model = anno.transcripts[tx_id]

            # For MXE, the first exon is 'inclusion' and second is 'skipped'.
            if self.ijc_sample_1 >= min_ijc:
                aln = first_downstream_junction.align_to_transcript(tx_model, True, False)
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)

            if self.sjc_sample_1 > min_sjc:
                aln = second_upstream_junction.align_to_transcript(tx_model, False, True)
                if aln:
                    variants += aln.convert_to_variant_records(anno, gene_seq, var_id)

        return list(set(variants))
