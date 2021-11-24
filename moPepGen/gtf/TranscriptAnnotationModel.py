""" Module for Transcript Annotation model """
from __future__ import annotations
from typing import List, TYPE_CHECKING
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen import dna, ERROR_INDEX_IN_INTRON


if TYPE_CHECKING:
    from Bio.Seq import Seq

GTF_FEATURE_TYPES = ['transcript', 'cds', 'exon', 'start_codon', 'stop_codon',
    'utr']

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
    def __init__(self, transcript:GTFSeqFeature=None,
            cds:List[GTFSeqFeature]=None, exon:List[GTFSeqFeature]=None,
            start_codon:List[GTFSeqFeature]=None,
            stop_codon:List[GTFSeqFeature]=None,
            utr:List[GTFSeqFeature]=None,
            five_utr:List[GTFSeqFeature]=None,
            three_utr:List[GTFSeqFeature]=None,
            is_protein_coding:bool=None):
        """ Construct a TranscriptAnnotationmodel """
        self.transcript = transcript
        self.cds = cds or []
        self.exon = exon or []
        self.start_codon = start_codon or []
        self.stop_codon = stop_codon or []
        self.utr = utr or []
        self.five_utr = five_utr or []
        self.three_utr = three_utr or []
        self.is_protein_coding = is_protein_coding

    def add_record(self, _type:str, record: GTFSeqFeature):
        """ Add a GTFRecrod into a TranscriptAnnotationModel. If the biotype
        is cds, exon, start_codon, or stop_codon, it is inserted in order.

        Args:
            type (str): Type of annotation to add. Must be from transcript,
                cds, exon, start_codon, or stop_codon.
            record (GTFRecord): The GTF record to be added.
        """
        if _type not in GTF_FEATURE_TYPES:
            raise ValueError(f'Type must be from {GTF_FEATURE_TYPES}')
        if _type == 'transcript':
            self.transcript = record
        else:
            if self.__getattribute__(_type) is None:
                self.__setattr__(_type, [])
            self.__getattribute__(_type).append(record)

    def split_utr(self) -> None:
        """ Infer 3'UTR or 5'UTR. CDS must be already sorted. """
        if not self.utr:
            return
        if not self.cds:
            raise ValueError(
                "UTR found but not CDS for transcript" + \
                    self.transcript.transcript_id
            )
        strand = self.transcript.strand
        for utr in self.utr:
            if strand == 1 and utr < self.cds[0] \
                    or strand == -1 and utr > self.cds[-1]:
                self.five_utr.append(utr)
            else:
                self.three_utr.append(utr)

    def sort_records(self):
        """ sort records """
        self.cds.sort()
        self.exon.sort()
        self.start_codon.sort()
        self.stop_codon.sort()
        self.split_utr()
        self.utr.sort()
        self.three_utr.sort()
        self.five_utr.sort()

    def get_cds_start_index(self) -> int:
        """ Returns the CDS start index of the transcript """
        cds_start = 0
        if self.transcript.strand == 1:
            cds_frame = self.cds[0].frame
            for exon in self.exon:
                if self.cds[0].location.start not in exon:
                    cds_start += len(exon)
                    continue
                cds_start += (self.cds[0].location.start - exon.location.start)
                break
        elif self.transcript.strand == -1:
            cds_frame = self.cds[-1].frame or 0
            for exon in reversed(self.exon):
                if self.cds[-1].location.end != exon.location.end and \
                        self.cds[-1].location.end not in exon:
                    cds_start += len(exon)
                    continue
                cds_start += (exon.location.end - self.cds[-1].location.end)
                break
        else:
            raise ValueError('Strand must not be unknown.')
        return cds_start + cds_frame

    def get_cds_end_index(self, seq:Seq, start:int) -> int:
        """ Returns the CDS stop index of the transcript. """
        if self.three_utr:
            end = self.get_transcript_index(self.three_utr[0].location.start)
            return end - (end - start) % 3
        return len(seq) - (len(seq) - start) % 3

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

        if self.transcript.strand == -1:
            seq = seq.reverse_complement()

        if self.cds:
            cds_start = self.get_cds_start_index()
            cds_end = self.get_cds_end_index(seq, cds_start)
            orf = FeatureLocation(start=cds_start, end=cds_end)
        else:
            orf = None

        location = MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.transcript_id,
                start=0,
                end=len(seq)
            )
        )
        transcript = dna.DNASeqRecordWithCoordinates(
            seq=seq,
            orf=orf,
            locations=[location]
        )
        transcript.id = self.transcript.transcript_id
        transcript.name = self.transcript.transcript_id
        transcript.description = self.transcript.transcript_id +\
            '|' + self.transcript.gene_id
        if 'protein_id' in self.transcript.attributes:
            transcript.description += '|' + \
                self.transcript.protein_id
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

        location = MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.transcript_id,
                start=cds_start,
                end=cds_start + len(seq)
            )
        )
        cdna = dna.DNASeqRecordWithCoordinates(
            seq=seq,
            locations=[location]
        )
        cdna.id = self.transcript.transcript_id
        cdna.name = self.transcript.transcript_id
        cdna.description = self.transcript.transcript_id + '|' \
            + self.transcript.gene_id + '|' \
            + self.transcript.protein_id
        return cdna

    def get_transcript_index(self, genomic_index:int) -> int:
        """ Get the transcript index from a genomic index. """
        if self.transcript.strand == 1:
            index = 0
            if genomic_index < self.exon[0].location.start \
                    or genomic_index > self.exon[-1].location.end:
                raise ValueError(
                    r'The genomic index isn\'t in the range of this transcript'
                )
            for exon in self.exon:
                if exon.location.end < genomic_index:
                    index += exon.location.end - exon.location.start
                elif exon.location.start <= genomic_index:
                    index += genomic_index - exon.location.start
                    break
                else:
                    raise ValueError(ERROR_INDEX_IN_INTRON)
        elif self.transcript.strand == -1:
            index = -1
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

    def is_cds_start_nf(self) -> bool:
        """ Returns if the transcript has the tag of cds_start_NF """
        return 'tag' in self.transcript.attributes and \
            'cds_start_NF' in self.transcript.attributes['tag']

    def transcript_len(self) -> int:
        """ Get the transcript length minus introns """
        length = 0
        for exon in self.exon:
            length += exon.location.end - exon.location.start
        return length

    def say_hello(self):
        print('hello')
