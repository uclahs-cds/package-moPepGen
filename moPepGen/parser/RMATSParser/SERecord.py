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
            downstream_exon_start:int, downstream_exon_end:int):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.upstream_exon_start = upstream_exon_start
        self.upstream_exon_end = upstream_exon_end
        self.downstream_exon_start = downstream_exon_start
        self.downstream_exon_end = downstream_exon_end

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
            downstream_exon_end=int(fields[10])
        )

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
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
        end = anno.coordinate_genomic_to_gene(self.exon_end, self.gene_id)
        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        ref = str(gene_seq.seq[start])

        if anno.genes[self.gene_id].location.strand == -1:
            start, end = end, start

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        if not skipped:
            for tx_id in retained:
                alt = '<DEL>'
                attrs = {
                    'TRANSCRIPTS': tx_id,
                    'START': start,
                    'END': end,
                    'GENE_SYMBOL': gene_model.attributes['gene_name'],
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Deletion'
                _id = f'SE_{start}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not retained:
            if gene_model.location.strand == 1:
                insert_position = anno.coordinate_genomic_to_gene(
                    self.upstream_exon_end, self.gene_id
                )
            else:
                insert_position = anno.coordinate_genomic_to_gene(
                    self.downstream_exon_start, self.gene_id
                )
            for tx_id in skipped:
                location = FeatureLocation(
                    seqname=tx_id,
                    start=insert_position,
                    end=insert_position + 1
                )
                ref = str(gene_seq.seq[insert_position])
                alt = '<INS>'
                attrs = {
                    'TRANSCRIPTS': tx_id,
                    'DONOR_START': start,
                    'DONOR_END': end,
                    'DONOR_GENE_ID': self.gene_id,
                    'GENE_SYMBOL': gene_model.attributes['gene_name'],
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Insertion'
                _id = f'SE_{start}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        return variants
