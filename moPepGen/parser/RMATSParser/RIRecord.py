""" Module for Retained Intron """
from __future__ import annotations
from typing import List, Tuple
from moPepGen import seqvar, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord


class RIRecord(RMATSRecord):
    """ Retained Intron """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            retained_intron_exon_start:int, retained_intron_exon_end:int,
            upstream_exon_start:int, upstream_exon_end:int,
            downstream_exon_start:int, downstream_exon_end:int):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.retained_intron_exon_start = retained_intron_exon_start
        self.retained_intron_exon_end = retained_intron_exon_end
        self.upstream_exon_start = upstream_exon_start
        self.upstream_exon_end = upstream_exon_end
        self.downstream_exon_start = downstream_exon_start
        self.downstream_exon_end = downstream_exon_end

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
            downstream_exon_end=int(fields[10])
        )

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        """ Convert to list of VariantRecord """
        variants = []
        transcript_ids = anno.genes[self.gene_id].transcripts
        chrom = anno.genes[self.gene_id].location.seqname

        have_adjacent:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []

        for transcript_id in transcript_ids:
            model = anno.transcripts[transcript_id]
            it = iter(model.exon)
            exon = next(it, None)

            while exon:
                if int(exon.location.start) == self.upstream_exon_start and \
                        int(exon.location.end) == self.upstream_exon_end:
                    exon = next(it, None)
                    if exon and \
                            exon.location.start == self.downstream_exon_start and \
                            exon.location.end == self.downstream_exon_end:
                        have_adjacent.append((transcript_id, model))
                        break
                exon = next(it, None)

        start_gene = anno.coordinate_genomic_to_gene(
            self.upstream_exon_end, self.gene_id)
        end_gene = anno.coordinate_genomic_to_gene(
            self.downstream_exon_start, self.gene_id)

        if anno.genes[self.gene_id].location.strand == -1:
            start_gene, end_gene = end_gene, start_gene

        genomic_position = f'{chrom}:{self.upstream_exon_end}-{self.downstream_exon_start}'

        for transcript_id, model in have_adjacent:
            position = anno.coordinate_gene_to_transcript(start_gene,
                self.gene_id, transcript_id)

            location = FeatureLocation(seqname=transcript_id, start=position,
                end=position + 1)
            seq = model.get_transcript_sequence(genome[chrom])
            ref = str(seq.seq[position])
            alt = '<INS>'
            attrs = {
                'GENE_ID': self.gene_id,
                'START': start_gene,
                'END': end_gene,
                'COORDINATE': 'gene',
                'GENE_SYMBOL': model.transcript.attributes['gene_name'],
                'GENOMIC_POSITION': genomic_position
            }
            _type = 'Insertion'
            _id = f'RI_{start_gene}'
            record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
            variants.append(record)
        return variants
