""" Model for DNA sequence data
"""
from typing import Union, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class DNASeqRecord(SeqRecord):
    """ A DNASeqRecord object holds a single DNA sequence and information about
    it. Derived from the Bio.SeqRecord.SeqRecord class.
    """

class DNASeqDict(dict):
    """ A DNASeqDict object is a dict-like object that the values are
    DNASeqRecord onbjects.
    """
    def __init__(self, *args, **kwargs)->None:
        for val in kwargs.values():
            self._validate(val)
        super().__init__(*args, **kwargs)
    
    @staticmethod
    def _validate(val):
        if not isinstance(val, DNASeqRecord):
            raise TypeError(
                "'DNASeqDict' only accepts 'DNASeqRecord' objects."
            )

    def __setitem__(self, k:str, v:DNASeqRecord)->None:
        """ Set items. Only allow DNASeqRecord in values. """
        self._validate(v)
        super().__setitem__(k, v)
    
    def dump_genome(self, path:str)->None:
        """ Dump a FASTA file to a DNASeqDict
        
        Args:
            path (str): Path to the genome assembly FASTA file.
        """
        for record in SeqIO.parse(path, 'fasta'):
            record.__class__ = DNASeqRecord
            if record.id in self.keys():
                raise ValueError(
                    'Duplicated seqnames found in FASTA file: ' + path
                )
            self[record.id] = record