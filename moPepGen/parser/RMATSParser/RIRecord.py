""" Module for Retained Intron """
from __future__ import annotations
from typing import List
from moPepGen import seqvar, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord


class RIRecord(RMATSRecord):
    """ Retained Intron """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            retained_intron_exon_start:int, retained_intron_exon_end:int,
            upstream_exon_start:int, upstream_exon_end:int,
            downstream_exon_start:int, downstream_exon_end:int,
            ijc_sample_1:int, sjc_sample_1:int, ijc_sample_2:int,
            sjc_sample_2:int, inc_form_len:int, skip_form_len:int,
            pvalue:float, fdr:float):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.retained_intron_exon_start = retained_intron_exon_start
        self.retained_intron_exon_end = retained_intron_exon_end
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
    def readline(cls, line:str) -> RIRecord:
        """ Create a A5SSRecord from a line of the rMASTS output """
        fields = line.rstrip().split('\t')
        return cls(
            gene_id=fields[1].strip('"'),
            gene_symbol=fields[2].strip('"'),
            chrom=fields[3],
            retained_intron_exon_start=int(fields[5]),
            retained_intron_exon_end=int(fields[6]),
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

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        """ Convert to list of VariantRecord """
        variants = []
        gene_model = anno.genes[self.gene_id]
        transcript_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])

        spliced_in_ref:List[str] = []
        retained_in_ref:List[str] = []

        for tx_id in transcript_ids:
            model = anno.transcripts[tx_id]
            it = iter(model.exon)
            exon = next(it, None)

            while exon:
                if int(exon.location.end) == self.upstream_exon_end:
                    exon = next(it, None)
                    if not exon:
                        break
                    if int(exon.location.start) == self.downstream_exon_start:
                        spliced_in_ref.append(tx_id)

                        break
                exon_start = int(exon.location.start)
                exon_end = int(exon.location.end)

                if exon_start < self.upstream_exon_end \
                        < self.downstream_exon_start < exon_end - 1:
                    retained_in_ref.append(tx_id)
                exon = next(it, None)

        start_gene = anno.coordinate_genomic_to_gene(
            self.upstream_exon_end, self.gene_id)
        end_gene = anno.coordinate_genomic_to_gene(
            self.downstream_exon_start - 1, self.gene_id)

        if gene_model.location.strand == -1:
            start_gene, end_gene = end_gene, start_gene
        end_gene += 1

        genomic_position = f'{chrom}:{self.upstream_exon_end}-{self.downstream_exon_start}'

        var_id = f"RI_{start_gene}-{end_gene}"

        if not retained_in_ref and self.ijc_sample_1 >= min_ijc:
            insert_position = start_gene - 1
            location = FeatureLocation(seqname=self.gene_id,
                start=insert_position, end=insert_position + 1)
            for tx_id in spliced_in_ref:
                ref = str(gene_seq.seq[insert_position])
                alt = '<INS>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'DONOR_START': start_gene,
                    'DONOR_END': end_gene,
                    'DONOR_GENE_ID': self.gene_id,
                    'COORDINATE': 'gene',
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Insertion'
                _id = var_id
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        if not spliced_in_ref and self.sjc_sample_1 >= min_sjc:
            del_start = start_gene
            del_end = end_gene
            location = FeatureLocation(seqname=self.gene_id,
                start=del_start, end=del_end)
            for tx_id in retained_in_ref:
                ref = str(gene_seq.seq[del_start])
                alt = '<DEL>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'START': del_start,
                    'END': del_end,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Deletion'
                _id = var_id
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        return variants
