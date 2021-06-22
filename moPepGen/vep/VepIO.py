""" The IO module for VEP output files
"""
from typing import IO, Union, Iterable
from Bio.SeqIO.Interfaces import SequenceIterator
from moPepGen import vep


class VepIterator(SequenceIterator):
    """ Iterator for reading a VEP file.
    """
    def __init__(self, source:Union[IO, str], mode='t'):
        """ Constructor """
        super().__init__(source=source, mode=mode, fmt='VEP')

    def parse(self, handle:IO[str]) -> Iterable[vep.VEPRecord]:
        """ parse """
        records = self.iterate(handle)
        return records

    @staticmethod
    def iterate(handle:IO[str]) -> vep.VEPRecord:
        """ Read the VEP output tsv file and returns a generator """
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            fields = line.split('\t')

            consequences = fields[6].split(',')

            amino_acids = tuple(aa for aa in fields[10].split('/'))
            if len(amino_acids) == 1:
                amino_acids = (amino_acids[0], '')

            codons = tuple(codon for codon in fields[11].split('/'))
            if len(codons) == 1:
                codons = (codons[0], '')

            extra = {}
            for field in fields[13].split(';'):
                key, val = field.split('=')
                extra[key] = val

            yield vep.VEPRecord(
                uploaded_variation=fields[0],
                location=fields[1],
                allele=fields[2],
                gene=fields[3],
                feature=fields[4],
                feature_type=fields[5],
                consequences=consequences,
                cdna_position=fields[7],
                cds_position=fields[8],
                protein_position=fields[9],
                amino_acids=amino_acids,
                codons=codons,
                existing_variation='' if fields[12] == '-' else fields[12],
                extra=extra
            )


def parse(handle:Union[IO[str], str]) -> VepIterator:
    """ parse a VEP tsv file """
    return VepIterator(handle)

def read(handle:Union[IO[str], str]) -> vep.VEPRecord:
    """ read a vep tsv file """
    iterator = parse(handle)
    try:
        record = next(iterator)
    except StopIteration as e:
        raise ValueError('No records found in handle') from e

    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    return record
