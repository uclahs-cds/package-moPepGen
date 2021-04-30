""" This module defines the class logic for the VEP variants.
"""
from typing import List, Iterable
from moPepGen.io import VepIO
from moPepGen.VEPRecord import VEPRecord


class TranscriptVEPRecords(list):
    """ A TranscriptVEPRecords holds all VEP records associated with the same
    transcript ID.


    """
    def __init__(self, *args):
        for it in args:
            if not isinstance(val, VEPRecord):
                raise ValueError('Items must be VEPRecord')
        super().__init__(*args)

    @staticmethod
    def _validate(val):
        """ Validate parameters """
        if not isinstance(val, VEPRecord):
            raise TypeError(
                '\'TranscriptVEPRecords\' object only accepts \'VEPRecord\''
                'objects.'
            )

    def __setitem__(self, key:List[int], val:Iterable[VEPRecord]):
        """ Set item """
        for v in val:
            self._validate(v)
        super().__setitem__(key, val)

    def append(self, object:VEPRecord):
        """ Add a VEPRecord item to the end of the list. """
        self._validate(object)
        super().append(object)
    
    def extend(self, iterable:Iterable[VEPRecord]):
        """ Extend VEPRecords """
        for val in iterable:
            self._validate(val)
        super().extend(iterable)
    
    def insert(self, index:int, object:VEPRecord):
        """ Insert a VEPRecord at given position """
        self._validate(object)
        super().insert(index, object)

class TranscriptVEPDict(dict):
    """ A TranscriptGTFDict object is a dict-like object that holds the VEP
    results, that the keys are the transcript IDs and values are instances of
    VEPRecord.
    """
    def __init__(self, *args, **kwargs):
        for val in kwargs.values():
            if not isinstance(val, TranscriptVEPRecords):
                raise ValueError(
                    'The values of a TranscriptGTFDict must be '
                    'TranscriptAnnotationModel.'
                )
        super().__init__(*args, **kwargs)
    
    @staticmethod
    def _validate(val):
        """ Validate parameters """
        if not isinstance(val, TranscriptVEPRecords):
            raise TypeError(
                "'TranscriptVEPDict' object only accepts "
                'TranscriptVEPRecords objects.'
            )

    def __setitem__(self, k: str, v: TranscriptVEPRecords):
        """ set item """
        self._validate(v)
        super().__setitem__(k, v)

    def dump_vep(self, path:str)->None:
        """ Dump a VEP file to a TranscriptVEPDict
        
        Args:
            path (str): Path to the VEP output file.
        """
        for record in VepIO.parse(path):
            transcript_id = record.feature
            if transcript_id not in self.keys():
                self[transcript_id] = TranscriptVEPRecords()
            self[transcript_id].append(record)
    