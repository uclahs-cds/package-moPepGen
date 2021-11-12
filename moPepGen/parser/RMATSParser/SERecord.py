""" Module for rMATS parser """
from __future__ import annotations
from typing import List
from moPepGen import gtf, dna, seqvar
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord

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

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        """ Convert to variant records """
        variants = []
        gene_model = anno.genes[self.gene_id]
        transcript_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])

        skipped:List[str] = []
        retained:List[str] = []

        for tx_id in transcript_ids:
            model = anno.transcripts[tx_id]
            it = iter(model.exon)
            exon = next(it, None)

            while exon:
                if exon.location.start != self.upstream_exon_start:
                    exon = next(it, None)
                    continue

                exon = next(it, None)
                if not exon:
                    continue

                if int(exon.location.start) == self.downstream_exon_start:
                    skipped.append(tx_id)
                    break
                if exon.location.start == self.exon_start and \
                        exon.location.end == self.exon_end:
                    exon = next(it, None)
                    if not exon:
                        continue
                    if int(exon.location.start) == self.downstream_exon_start:
                        retained.append(tx_id)
                        break

        if skipped and retained:
            return variants

        if not skipped and not retained:
            return variants

        start = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        end = anno.coordinate_genomic_to_gene(self.exon_end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        # For SE, the skipped is 'skipped', and inclusion is 'inclusion'
        # what a nonsense comment
        if not skipped and self.sjc_sample_1 >= min_sjc:
            location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
            for tx_id in retained:
                alt = '<DEL>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'START': start,
                    'END': end,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Deletion'
                _id = f'SE_{start}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not retained and self.ijc_sample_1 >= min_ijc:
            if gene_model.strand == 1:
                insert_position = anno.coordinate_genomic_to_gene(
                    index=self.upstream_exon_end - 1,
                    gene=self.gene_id
                )
            else:
                insert_position = anno.coordinate_genomic_to_gene(
                    index=self.downstream_exon_start,
                    gene=self.gene_id
                )
            location = FeatureLocation(seqname=self.gene_id, start=insert_position,
                end=insert_position + 1)
            for tx_id in skipped:
                ref = str(gene_seq.seq[insert_position])
                alt = '<INS>'
                attrs = {
                    'TRANSCRIPT_ID': tx_id,
                    'DONOR_START': start,
                    'DONOR_END': end,
                    'DONOR_GENE_ID': self.gene_id,
                    'GENE_SYMBOL': gene_model.gene_name,
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Insertion'
                _id = f'SE_{start}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        return variants
