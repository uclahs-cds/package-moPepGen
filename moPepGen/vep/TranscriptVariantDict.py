""" This module defines the class logic for the VEP variants.
"""
from typing import List, Iterable
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord


class VEPVariantRecords(list):
    """ A TranscriptVEPRecords holds all VEP records associated 
    with the same transcript ID.
    """
    def __init__(self, *args):
        for it in args:
            if not isinstance(it, VEPVariantRecord):
                raise ValueError('Items must be VEPVariantRecord')
        super().__init__(*args)

    @staticmethod
    def _validate(val):
        """ Validate parameters """
        if not isinstance(val, VEPVariantRecord):
            raise TypeError(
                '\'VEPVariantRecord\' object only accepts' 
                '\'VEPVariantRecord\' objects.'
            )

    def __setitem__(self, key:List[int], val:Iterable[VEPVariantRecord]):
        """ Set item """
        for v in val:
            self._validate(v)
        super().__setitem__(key, val)

    def append(self, object:VEPVariantRecord):
        """ Add a VEPRecord item to the end of the list. """
        self._validate(object)
        super().append(object)
    
    def extend(self, iterable:Iterable[VEPVariantRecord]):
        """ Extend VEPRecords """
        for val in iterable:
            self._validate(val)
        super().extend(iterable)
    
    def insert(self, index:int, object:VEPVariantRecord):
        """ Insert a VEPVariantRecord at given position """
        self._validate(object)
        super().insert(index, object)


class TranscriptVariantDict(dict):
    """ A TranscriptVariantDict object is a dict-like object that holds the
    variant results for each transcipts. The keys are the transcript IDs and
    values are list of VEPVariantRecord.
    """
    
    def _validate(self, val):
        """ Validate parameters """
        if not isinstance(val, VEPVariantRecords):
            raise TypeError(
                "'TranscriptVEPDict' object only accepts "
                'TranscriptVEPRecords objects.'
            )

    def __setitem__(self, k: str, v: List[VEPVariantRecord]):
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
                result += f"...\n"
                i = len(self) - 4
            i += 1
        result += f'\n{len(self)} transcripts'
        return result
    