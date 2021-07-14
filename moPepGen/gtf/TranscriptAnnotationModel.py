""" Module for Transcript Annotation model """
from typing import List
from moPepGen.SeqFeature import SeqFeature, FeatureLocation
from moPepGen import dna


_GTF_FEATURE_TYPES = ['transcript', 'cds', 'exon', 'start_codon', 'stop_codon']
INDEX_IN_INTRON_ERROR = 'The genomic index seems to be in an intron'

class TranscriptAnnotationModel():
    """ A TranscriptAnnotationModel holds all the annotations associated with
    the same transcript from a GTF file.

    Attributes:
        transcript (GTFRecord): The GTF record of the transcript
        cds (List[GTFRecord]): The GTF records of all the coding sequences of
            the transcript.
        exon (List[GTFRecord]): The GTF records of all the exons of the
            transcript.
        start_codon (List[GTFRecord]): The GTF records of all the start codons
            of the transcript.
        stop_codon (List[GTFRecord]): The GTF records of all the stop codons
            of the transcript.
    """
    def __init__(self, transcript:SeqFeature=None,
            cds:List[SeqFeature]=None, exon:List[SeqFeature]=None,
            start_codon:List[SeqFeature]=None,
            stop_codon:List[SeqFeature]=None):
        """ Construct a TranscriptAnnotationmodel """
        self.transcript = transcript
        self.cds = [] if cds is None else cds
        self.exon = [] if exon is None else exon
        self.start_codon = [] if start_codon is None else start_codon
        self.stop_codon = [] if stop_codon is None else stop_codon

    def add_record(self, _type:str, record: SeqFeature):
        """ Add a GTFRecrod into a TranscriptAnnotationModel. If the biotype
        is cds, exon, start_codon, or stop_codon, it is inserted in order.

        Args:
            type (str): Type of annotation to add. Must be from transcript,
                cds, exon, start_codon, or stop_codon.
            record (GTFRecord): The GTF record to be added.
        """
        if _type not in _GTF_FEATURE_TYPES:
            raise ValueError(f'Type must be from {_GTF_FEATURE_TYPES}')
        if _type == 'transcript':
            self.transcript = record
        else:
            if self.__getattribute__(_type) is None:
                self.__setattr__(_type, [])
            self.__getattribute__(_type).append(record)

    def sort_records(self):
        """ sort records """
        self.cds.sort()
        self.exon.sort()
        self.start_codon.sort()
        self.stop_codon.sort()

    def get_cds_start_index(self) -> int:
        """ Returns the CDS start index of the transcript """
        cds_start = 0
        if self.transcript.strand == 1:
            for exon in self.exon:
                if self.cds[0].location.start not in exon:
                    cds_start += len(exon)
                    continue
                cds_start += (self.cds[0].location.start - exon.location.start)
                break
        elif self.transcript.strand == -1:
            for exon in reversed(self.exon):
                if self.cds[-1].location.end != exon.location.end and \
                        self.cds[-1].location.end not in exon:
                    cds_start += len(exon)
                    continue
                cds_start += (exon.location.end - self.cds[-1].location.end)
                break
        else:
            raise ValueError('Strand must not be unknown.')
        return cds_start

    def get_cds_end_index(self) -> int:
        """ Returns the CDS stop index of the transcript. """
        cds_end = 0
        if self.transcript.strand == 1:
            for exon in self.exon:
                if self.cds[-1].location.end not in exon:
                    cds_end += len(exon)
                    continue
                cds_end += self.cds[-1].location.end - exon.location.start
                break
        elif self.transcript.strand == -1:
            for exon in reversed(self.exon):
                if self.cds[0].location.start not in exon:
                    cds_end += len(exon)
                    continue
                cds_end += exon.location.end - self.cds[0].location.start
        else:
            raise ValueError('Strand must not be unknown.')
        return cds_end

    def get_transcript_sequence(self, chrom:dna.DNASeqRecord
            ) -> dna.DNASeqRecordWithCoordinates:
        """ Returns the transcript sequence. The is done by concating all the
        exon sequences. If the gene is on the negatice strand, the reverse
        complement is returned.

        Args:
            chrom (DNASeqRecord): The chromosome sequence that the transcript
                is located.
        """
        if len(self.exon) == 0:
            raise ValueError("Transcript model has no exon")
        seq = None
        for location in [it.location for it in self.exon]:
            new_seq = chrom.seq[location.start:location.end]
            if seq is None:
                seq = new_seq
            else:
                seq = seq + new_seq

        if self.cds:
            cds_start = self.get_cds_start_index()
            cds_end = self.get_cds_end_index()
            orf = FeatureLocation(start=cds_start, end=cds_end)
        else:
            orf = None

        if self.transcript.strand == -1:
            seq = seq.reverse_complement()

        location = dna.MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.attributes['transcript_id'],
                start=0,
                end=len(seq)
            )
        )
        transcript = dna.DNASeqRecordWithCoordinates(
            seq=seq,
            orf=orf,
            locations=[location]
        )
        transcript.id = self.transcript.attributes['transcript_id']
        transcript.name = self.transcript.attributes['transcript_id']
        transcript.description = self.transcript.attributes['transcript_id'] +\
            '|' + self.transcript.attributes['gene_id']
        if 'protein_id' in self.transcript.attributes:
            transcript.description += '|' + \
                self.transcript.attributes['protein_id']
        return transcript

    def get_cdna_sequence(self, chrom:dna.DNASeqRecord
            ) -> dna.DNASeqRecordWithCoordinates:
        """ Get the cDNA sequence.

        Args:
            chrom (DNASeqRecord): The chromosome sequence.

        Returns:
            The cDNA sequence (DNASeqRecordWithCoordinates)
        """
        if len(self.cds) == 0:
            raise ValueError('Transcript model has no cds')
        seq = None
        locations = [it.location for it in self.cds]
        for location in locations:
            new_seq = chrom.seq[location.start:location.end]
            if seq is None:
                seq = new_seq
            else:
                seq = seq + new_seq

        cds_start = self.get_cds_start_index()

        if self.transcript.strand == -1:
            seq = seq.reverse_complement()

        location = dna.MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.attributes['transcript_id'],
                start=cds_start,
                end=cds_start + len(seq)
            )
        )
        cdna = dna.DNASeqRecordWithCoordinates(
            seq=seq,
            locations=[location]
        )
        cdna.id = self.transcript.attributes['transcript_id']
        cdna.name = self.transcript.attributes['transcript_id']
        cdna.description = self.transcript.attributes['transcript_id'] + '|' \
            + self.transcript.attributes['gene_id'] + '|' \
            + self.transcript.attributes['protein_id']
        return cdna

    def get_transcript_index(self, genomic_index:int) -> int:
        """ Get the transcript index from a genomic index. """
        index = 0
        if self.transcript.strand == 1:
            if genomic_index < self.exon[0].location.start \
                    or genomic_index > self.exon[-1].location.end:
                raise ValueError(
                    r'The genomic index isn\'t in the range of this transcript'
                )
            for exon in self.exon:
                if exon.location.end < genomic_index:
                    index += exon.location.end - exon.location.start
                elif exon.location.start < genomic_index:
                    index += genomic_index - exon.location.start
                    break
                else:
                    raise ValueError(INDEX_IN_INTRON_ERROR)
        elif self.transcript.strand == -1:
            if genomic_index < self.exon[0].location.start \
                    or genomic_index > self.exon[-1].location.end:
                raise ValueError(
                    r'The genomic index isn\'t in the range of this transcript'
                )
            for exon in reversed(self.exon):
                if exon.location.start > genomic_index:
                    index += exon.location.end - exon.location.start
                elif exon.location.end > genomic_index:
                    index += exon.location.end - genomic_index
                    break
                else:
                    raise ValueError(
                        'The genomic index seems to be in an intron'
                    )
        return index

    def get_transcript_start_genomic_coordinate(self) -> int:
        """ Get the genomic coordinate of the start point of the transcript """
        if self.transcript.location.strand == 1:
            return int(self.exon[0].location.start)
        return int(self.exon[-1].location.end)
