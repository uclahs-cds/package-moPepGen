""" Model for DNA sequence data
"""
from Bio import SeqIO
from moPepGen.dna import DNASeqRecord


class DNASeqDict(dict):
    """ A DNASeqDict object is a dict-like object that the values are
    DNASeqRecord objects.
    """
    def __init__(self, *args, **kwargs) -> None:
        for val in kwargs.values():
            self._validate(val)
        super().__init__(*args, **kwargs)

    @staticmethod
    def _validate(val):
        """ validate values """
        if not isinstance(val, DNASeqRecord):
            raise TypeError(
                "'DNASeqDict' only accepts 'DNASeqRecord' objects."
            )

    def __setitem__(self, k:str, v:DNASeqRecord) -> None:
        """ Set items. Only allow DNASeqRecord in values. """
        self._validate(v)
        super().__setitem__(k, v)

    def dump_fasta(self, path:str) -> None:
        """ Dump a FASTA file to a DNASeqDict

        Args:
            path (str): Path to the genome assembly FASTA file.
        """
        for record in SeqIO.parse(path, 'fasta'):
            record.__class__ = DNASeqRecord
            record.id = record.id.split(' ')[0]
            record.id = record.id.split('|')[0]
            if record.id in self.keys():
                raise ValueError(
                    'Duplicated seqnames found in FASTA file: ' + path
                )
            self[record.id] = record

    def with_chr(self) -> bool:
        """ Check whether keys have 'chr' or not """
        # so far only works with humans
        return bool('chr1' in self.keys())
