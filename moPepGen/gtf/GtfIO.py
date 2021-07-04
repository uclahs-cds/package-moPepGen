""" Module for GTF IO """
from typing import IO, Union, Iterable
from Bio.SeqIO.Interfaces import SequenceIterator
from moPepGen.SeqFeature import SeqFeature, FeatureLocation


class GtfIterator(SequenceIterator):
    """ GTF Iterator """
    def __init__(self, source:Union[IO, str], mode='t'):
        """ Constructor """
        super().__init__(source=source, mode=mode, fmt='GTF')

    def parse(self, handle:IO[str]) -> Iterable[SeqFeature]:
        """ parse
        """
        records = self.iterate(handle)
        return records

    @staticmethod
    def iterate(handle:IO[str]) -> Iterable[SeqFeature]:
        """ Iterate through a GTF file and yield a record each time. """
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            fields = line.split('\t')

            try:
                strand={'+':1, '-':-1, '?':0}[fields[6]]
            except KeyError:
                strand=None

            attributes = {}
            attribute_list = [field.strip().split(' ') for field in \
                fields[8].rstrip(';').split(';')]
            for key,val in attribute_list:
                val = val.strip('"')
                if key == 'tag':
                    if key not in attributes:
                        attributes[key] = []
                    attributes[key].append(val)
                else:
                    attributes[key] = val

            location=FeatureLocation(
                seqname=fields[0],
                start=int(fields[3])-1,
                end=int(fields[4])
            )

            record = SeqFeature(
                chrom=fields[0],
                attributes=attributes,
                location=location,
                type=fields[2],
                strand=strand
            )
            # record.id = record.attributes['transcript_id']
            yield record

def parse(handle:Union[IO[str], str]) -> GtfIterator:
    """ Parser for GTF files.
    """
    return GtfIterator(handle)
