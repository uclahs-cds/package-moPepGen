""" This module defines the class logic for the GTF annotations.
"""
from typing import List
from moPepGen.io import GtfIO
from moPepGen.GTFRecord import GTFRecord


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
            of the transcript
        stop_codon (List[GTFRecord]): The GTF records of all the stop codons
            of the transcript
    """
    def __init__(self, transcript:GTFRecord=None,
            cds:List[GTFRecord]=None, exon:List[GTFRecord]=None,
            start_codon:List[GTFRecord]=None, stop_codon:List[GTFRecord]=None):
        """ Construct a TranscriptAnnotationmodel """
        self.transcript = transcript
        self.cds = [] if cds is None else cds
        self.exon = [] if exon is None else exon
        self.start_codon = [] if start_codon is None else start_codon
        self.stop_codon = [] if stop_codon is None else stop_codon
    
    def add_record(self, type:str, record: GTFRecord):
        """ Add a GTFRecrod into a TranscriptAnnotationModel.

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
            raise TypeValue(
                'The values of a TranscriptGTFDict must be '
                'TranscriptAnnotationModel.'
            )

    def __setitem__(self, k:str, v:TranscriptAnnotationModel)->None:
        """ set item """
        self._validate(v)
        super().__setitem__(k, v)
    
    def dump_gtf(self, path:str)->None:
        """ Dump a GTF file into a TranscriptGTFDict
        
        Args:
            path (str): Path to a GTF file.
        """
        for record in GtfIO.parse(path):
            feature = record.feature.lower()
            if feature not in _GTF_FEATURE_TYPES:
                continue
            transcript_id = record.attributes['transcript_id']
            if transcript_id not in self.keys():
                self[transcript_id] = TranscriptAnnotationModel()
            self[transcript_id].add_record(feature, record)