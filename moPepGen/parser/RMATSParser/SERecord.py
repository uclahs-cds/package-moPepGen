""" Module for rMATS parser """
from __future__ import annotations
from typing import List, Tuple
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
        transcript_ids = anno.genes[self.gene_id].transcripts
        chrom = anno.genes[self.gene_id].location.seqname

        skipped:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []
        retained:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []

        for transcript_id in transcript_ids:
            model = anno.transcripts[transcript_id]
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
                    skipped.append((transcript_id, model))
                    break
                if exon.location.start == self.exon_start and \
                        exon.location.end == self.exon_end:
                    exon = next(it, None)
                    if not exon:
                        continue
                    if int(exon.location.start) == self.downstream_exon_start:
                        retained.append((transcript_id, model))
                        break

        if skipped and retained:
            return variants

        if not skipped and not retained:
            return variants

        start_gene = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        end_gene = anno.coordinate_genomic_to_gene(self.exon_end, self.gene_id)

        if anno.genes[self.gene_id].location.strand == -1:
            start_gene, end_gene = end_gene, start_gene

        genomic_position = f'{chrom}:{self.exon_start}:{self.exon_end}'\
            if model.transcript.location.strand == 1 else \
            f'{chrom}:{self.exon_end}-{self.exon_start}'

        if not skipped:
            for transcript_id, model in retained:
                start = anno.coordinate_gene_to_transcript(start_gene,
                    self.gene_id, transcript_id)
                end = anno.coordinate_gene_to_transcript(end_gene,
                    self.gene_id, transcript_id)
                location = FeatureLocation(
                    seqname=transcript_id,
                    start=start, end=end
                )
                seq = model.get_transcript_sequence(genome[chrom])
                ref = str(seq.seq[start])
                alt = '<DEL>'
                attrs = {
                    'GENE_ID': self.gene_id,
                    'START': start,
                    'END': end,
                    'GENE_SYMBOL': model.transcript.attributes['gene_name'],
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Deletion'
                _id = f'SE_{start_gene}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not retained:
            if model.transcript.location.strand == 1:
                insert_position_gene = anno.coordinate_genomic_to_gene(
                    self.upstream_exon_end, self.gene_id
                )
            else:
                insert_position_gene = anno.coordinate_genomic_to_gene(
                    self.downstream_exon_start, self.gene_id
                )
            for transcript_id, model in skipped:
                insert_position = anno.coordinate_gene_to_transcript(
                    insert_position_gene, self.gene_id, transcript_id
                )
                location = FeatureLocation(
                    seqname=transcript_id,
                    start=insert_position,
                    end=insert_position + 1
                )
                seq = model.get_transcript_sequence(genome[chrom])
                ref = str(seq.seq[insert_position])
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
                _id = f'SE_{start_gene}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        return variants
