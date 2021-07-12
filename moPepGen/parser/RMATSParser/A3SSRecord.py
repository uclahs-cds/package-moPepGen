""" Module fro Alternative 3' Splicing Site """
from __future__ import annotations
from typing import List, Tuple
from moPepGen import seqvar, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord


class A3SSRecord(RMATSRecord):
    """ Alternative 3' Splicing Site """
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
            flanking_exon_end=int(fields[10])
        )

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        variants = []
        transcript_ids = anno.genes[self.gene_id].transcripts
        chrom = anno.genes[self.gene_id].location.seqname

        short:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []
        long:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []

        for transcript_id in transcript_ids:
            model = anno.transcripts[transcript_id]
            if anno.genes[self.gene_id].location.strand == 1:
                it = iter(model.exon)
            else:
                it = reversed(model.exon)
            exon = next(it, None)

            while exon:
                if int(exon.location.start) == self.flanking_exon_start:
                    exon = next(it, None)
                    if not exon:
                        break
                    if exon.location.start == self.long_exon_start and \
                        int(exon.location.end) == self.long_exon_end:
                        long.append((transcript_id, model))
                    elif int(exon.location.start) == self.short_exon_start and \
                        int(exon.location.end) == self.short_exon_end:
                        short.append((transcript_id, model))
                    break
                exon = next(it, None)

        if (long and short) or (not long and not short):
            return variants

        if anno.genes[self.gene_id].location.strand == 1:
            start_gene = anno.coordinate_genomic_to_gene(self.long_exon_start,
                self.gene_id)
            end_gene = anno.coordinate_genomic_to_gene(self.short_exon_start,
                self.gene_id)
            genomic_position = f'{chrom}:{self.long_exon_start}:{self.short_exon_start}'
        else:
            start_gene = anno.coordinate_genomic_to_gene(self.long_exon_end,
                self.gene_id)
            end_gene = anno.coordinate_genomic_to_gene(self.short_exon_end,
                self.gene_id)
            genomic_position = f'{chrom}:{self.short_exon_end}-{self.long_exon_end}'

        if not short:
            for transcript_id, model in long:
                start = anno.coordinate_gene_to_transcript(start_gene,
                    self.gene_id, transcript_id)
                end = anno.coordinate_gene_to_transcript(end_gene,
                    self.gene_id, transcript_id)
                location = FeatureLocation(seqname=transcript_id, start=start,
                    end=end)
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
                _id = f'A5SS_{start_gene}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not long:
            if model.transcript.location.strand == 1:
                insert_position_gene = anno.coordinate_genomic_to_gene(
                    self.short_exon_end, self.gene_id
                )
            else:
                insert_position_gene = anno.coordinate_genomic_to_gene(
                    self.short_exon_start, self.gene_id
                )
            for transcript_id, model in short:
                insert_position = anno.coordinate_gene_to_transcript(
                    insert_position_gene, self.gene_id, transcript_id
                )
                start = anno.coordinate_gene_to_transcript(start_gene,
                    self.gene_id, transcript_id)
                end = anno.coordinate_gene_to_transcript(end_gene,
                    self.gene_id, transcript_id)
                location = FeatureLocation(
                    seqname=transcript_id,
                    start=insert_position,
                    end=insert_position+1
                )
                seq = model.get_transcript_sequence(genome[chrom])
                ref = str(seq.seq[start])
                alt = '<INS>'
                attrs = {
                    'GENE_ID': self.gene_id,
                    'START': start,
                    'END': end,
                    'COORDINATE': 'gene',
                    'GENE_SYMBOL': model.transcript.attributes['gene_name'],
                    'GENOMIC_POSITION': genomic_position
                }
                _type = 'Insertion'
                _id = f'SE_{start_gene}'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)
        return variants
