""" This module defines the class logic for the GTF annotations.
"""
from typing import List
from moPepGen.gtf import GtfIO
from moPepGen.SeqFeature import SeqFeature, FeatureLocation
from moPepGen.dna import DNASeqRecord, DNASeqRecordWithCoordinates, \
    DNASeqRecordWithCoordinates, MatchedLocation


_GTF_FEATURE_TYPES = ['transcript', 'cds', 'exon', 'start_codon', 'stop_codon']


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
    
    def add_record(self, type:str, record: SeqFeature):
        """ Add a GTFRecrod into a TranscriptAnnotationModel. If the biotype
        is cds, exon, start_codon, or stop_codon, it is inserted in order.

        Args:
            type (str): Type of annotation to add. Must be from transcript,
                cds, exon, start_codon, or stop_codon.
            record (GTFRecord): The GTF record to be added.
        """
        if type not in _GTF_FEATURE_TYPES:
            raise ValueError(f'Type must be from {_GTF_FEATURE_TYPES}')
        if type == 'transcript':
            self.transcript = record
        else:
            if self.__getattribute__(type) is None:
                self.__setattr__(type, [])
            self.__getattribute__(type).append(record)
    
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
                cds_end += len(self.cds[-1])
                break
        elif self.transcript.strand == -1:
            for exon in reversed(self.exon):
                if self.cds[0].location.start not in exon:
                    cds_end += len(exon)
                    continue
                cds_end += len(self.cds[0])
        else:
            raise ValueError('Strand must not be unknown.')
        return cds_end
        
    def get_transcript_sequence(self, chrom:DNASeqRecord
            )->DNASeqRecordWithCoordinates:
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
        cds_start = self.get_cds_start_index()
        cds_end = self.get_cds_end_index()
        if self.transcript.strand == -1:
            seq = seq.reverse_complement()
        location = MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.attributes['transcript_id'],
                start=0,
                end=len(seq)
            )
        )
        orf = FeatureLocation(start=cds_start, end=cds_end)
        transcript = DNASeqRecordWithCoordinates(
            seq=seq,
            orf=[orf],
            locations=[location]
        )
        transcript.id = self.transcript.attributes['transcript_id']
        transcript.name = self.transcript.attributes['transcript_id']
        transcript.description = self.transcript.attributes['transcript_id'] + '|' \
            + self.transcript.attributes['gene_id'] + '|' \
            + self.transcript.attributes['protein_id']
        return transcript

    def get_cdna_sequence(self, chrom:DNASeqRecord)->DNASeqRecordWithCoordinates:
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
                seqname=self.transcript.attributes['transcript_id'],
                start=cds_start,
                end=cds_start + len(seq)
            )
        )
        cdna = DNASeqRecordWithCoordinates(
            seq=seq,
            locations=[location]
        )
        cdna.id = self.transcript.attributes['transcript_id']
        cdna.name = self.transcript.attributes['transcript_id']
        cdna.description = self.transcript.attributes['transcript_id'] + '|' \
            + self.transcript.attributes['gene_id'] + '|' \
            + self.transcript.attributes['protein_id']
        return cdna

class TranscriptGTFDict(dict):
    """ A VEPTranscripts object is a dict-like object that holds the GTF
    results, that the keys are the transcript IDs and values are instances of
    VEPRecord.
    """
    def __init__(self, *args, **kwargs):
        """ Construct a TranscriptGTFDict object """
        for val in kwargs.values():
            self._validate(val)
        super().__init__(*args, **kwargs)
    
    @staticmethod
    def _validate(val):
        if not isinstance(val, TranscriptAnnotationModel):
            raise TypeError(
                'The values of a TranscriptGTFDict must be '
                'TranscriptAnnotationModel.'
            )

    def __setitem__(self, k:str, v:TranscriptAnnotationModel)->None:
        """ set item """
        self._validate(v)
        super().__setitem__(k, v)
    
    def __repr__(self) -> str:
        """ Return a string representation """
        result = ""
        i = 0
        while i < len(self):
            key = list(self.keys())[i]
            result += f"'{key}': {self[key].__repr__()}\n"
            if i == 3 and len(self) > 7:
                result += "...\n"
                i = len(self) - 4
            i += 1
        result += f'\n{len(self)} transcripts'
        return result
    
    def dump_gtf(self, path:str, biotype:List[str]=None)->None:
        """ Dump a GTF file into a TranscriptGTFDict
        
        Args:
            path (str): Path to a GTF file.
            biotype (List[str]): The annotation biotype to keep. Features
                in a GTF can be annotated as protein_coding, miRNA, lncRNA,
                miRNA, etc.
        """
        for record in GtfIO.parse(path):
            if biotype is not None and record.biotype not in biotype:
                continue
            feature = record.type.lower()
            if feature not in _GTF_FEATURE_TYPES:
                continue
            transcript_id = record.attributes['transcript_id']
            record.id = transcript_id
            if transcript_id not in self.keys():
                self[transcript_id] = TranscriptAnnotationModel()
            self[transcript_id].add_record(feature, record)
        for transcript_model in self.values():
            transcript_model.sort_records()