""" Module fro Mutually Exclusive Exons """
from __future__ import annotations
from typing import List, Tuple
from moPepGen import seqvar, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord

class MXERecord(RMATSRecord):
    """ Mutually Exclusive Exons """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            first_exon_start:int, first_exon_end:int, second_exon_start:int,
            second_exon_end:int, upstream_exon_start:int,
            upstream_exon_end:int, downstream_exon_start:int,
            downstream_exon_end:int):
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
            downstream_exon_end=int(fields[12])
        )

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        variants = []
        transcript_ids = anno.genes[self.gene_id].transcripts
        chrom = anno.genes[self.gene_id].location.seqname
        strand = anno.genes[self.gene_id].location.strand

        have_both:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []
        have_first:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []
        have_second:List[Tuple[str, gtf.TranscriptAnnotationModel]] = []

        for transcript_id in transcript_ids:
            model = anno.transcripts[transcript_id]
            it = iter(model.exon)
            exon = next(it, None)
            while exon:
                if int(exon.location.end) != self.upstream_exon_end:
                    exon = next(it, None)
                    continue
                exon = next(it, None)
                if not exon:
                    break
                if int(exon.location.start) == self.first_exon_start and \
                        int(exon.location.end) == self.first_exon_end:
                    exon = next(it, None)
                    if not exon:
                        break
                    if int(exon.location.start) == self.downstream_exon_start:
                        if strand == 1:
                            have_first.append((transcript_id, model))
                        else:
                            have_second.append((transcript_id, model))
                    elif int(exon.location.start) == self.second_exon_start and \
                            int(exon.location.end) == self.second_exon_end:
                        exon = next(it, None)
                        if not exon:
                            break
                        if int(exon.location.start) == self.downstream_exon_start:
                            have_both.append((transcript_id, model))
                elif int(exon.location.start) == self.second_exon_start and \
                        int(exon.location.end) == self.second_exon_end:
                    exon = next(it, None)
                    if not exon:
                        break
                    if int(exon.location.start) == self.downstream_exon_start:
                        if strand == 1:
                            have_second.append((transcript_id, model))
                        else:
                            have_first.append((transcript_id, model))

        if (have_first and have_second) or \
                (not have_first and not have_second and not have_both):
            return variants

        if anno.genes[self.gene_id].location.strand == 1:
            first_start_gene = anno.coordinate_genomic_to_gene(
                self.first_exon_start, self.gene_id)
            first_end_gene = anno.coordinate_genomic_to_gene(
                self.first_exon_end, self.gene_id)
            second_start_gene = anno.coordinate_genomic_to_gene(
                self.second_exon_start, self.gene_id)
            second_end_gene = anno.coordinate_genomic_to_gene(
                self.second_exon_end, self.gene_id)
            first_genomic_position = f'{chrom}:{self.first_exon_start + 1}-{self.first_exon_end}'
            second_genomic_position = f'{chrom}:{self.second_exon_start + 1}-{self.second_exon_end}'
        else:
            first_start_gene = anno.coordinate_genomic_to_gene(
                self.second_exon_end, self.gene_id)
            first_end_gene = anno.coordinate_genomic_to_gene(
                self.second_exon_start, self.gene_id)
            second_start_gene = anno.coordinate_genomic_to_gene(
                self.first_exon_end, self.gene_id)
            second_end_gene = anno.coordinate_genomic_to_gene(
                self.first_exon_start, self.gene_id)
            second_genomic_position = f'{chrom}:{self.first_exon_start + 1}-{self.first_exon_end}'
            first_genomic_position = f'{chrom}:{self.second_exon_end + 1}-{self.second_exon_end}'

        _id = f'MXE_{first_start_gene + 1}-{first_end_gene}:' +\
            f'{second_start_gene}-{second_end_gene}'

        if not have_second:
            for transcript_id, model in have_first:
                first_start = anno.coordinate_gene_to_transcript(
                    first_start_gene, self.gene_id, transcript_id)
                first_end = anno.coordinate_gene_to_transcript(
                    first_end_gene, self.gene_id, transcript_id)
                location = FeatureLocation(seqname=transcript_id,
                     start=first_start, end=first_end)
                seq = model.get_transcript_sequence(genome[chrom])
                ref = str(seq.seq[first_start])
                alt = '<SUB>'
                attrs = {
                    'GENE_ID': self.gene_id,
                    'START': first_start,
                    'END': first_end,
                    'DONOR_START': second_start_gene,
                    'DONOR_END': second_end_gene,
                    'COORDINATE': 'gene',
                    'GENE_SYMBOL': model.transcript.attributes['gene_name'],
                    'GENOMIC_POSITION': f'{chrom}-{first_start_gene + 1}:{first_end_gene}-'
                    f'{second_start_gene + 1}:{second_end_gene}'
                }
                _type = 'Substitution'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

            for transcript_id, model in have_both:
                start = anno.coordinate_gene_to_transcript(
                    first_start_gene, self.gene_id, transcript_id)
                end = anno.coordinate_gene_to_transcript(
                    first_end_gene, self.gene_id, transcript_id)
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
                    'GENOMIC_POSITION': first_genomic_position
                }
                _type = 'Deletion'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        if not have_first:
            for transcript_id, model in have_second:
                second_start = anno.coordinate_gene_to_transcript(
                    second_start_gene, self.gene_id, transcript_id)
                second_end = anno.coordinate_gene_to_transcript(
                    second_end_gene, self.gene_id, transcript_id)
                location = FeatureLocation(seqname=transcript_id,
                     start=second_start, end=second_end)
                seq = model.get_transcript_sequence(genome[chrom])
                ref = str(seq.seq[second_start])
                alt = '<SUB>'
                attrs = {
                    'GENE_ID': self.gene_id,
                    'START': second_start,
                    'END': second_end,
                    'DONOR_START': first_start_gene,
                    'DONOR_END': first_end_gene,
                    'COORDINATE': 'gene',
                    'GENE_SYMBOL': model.transcript.attributes['gene_name'],
                    'GENOMIC_POSITION': f'{chrom}-{second_start_gene + 1}:{second_end_gene}-'
                    f'{first_start_gene + 1}:{first_end_gene}'
                }
                _type = 'Substitution'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

            for transcript_id, model in have_both:
                start = anno.coordinate_gene_to_transcript(
                    second_start_gene, self.gene_id, transcript_id)
                end = anno.coordinate_gene_to_transcript(
                    second_end_gene, self.gene_id, transcript_id)
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
                    'GENOMIC_POSITION': second_genomic_position
                }
                _type = 'Deletion'
                record = seqvar.VariantRecord(location, ref, alt, _type, _id, attrs)
                variants.append(record)

        return variants
