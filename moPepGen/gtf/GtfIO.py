from typing import IO, Union
from Bio.SeqIO.Interfaces import SequenceIterator
from moPepGen.SeqFeature import SeqFeature, FeatureLocation


class GtfIterator(SequenceIterator):
    """"""
    def __init__(self, source:Union[IO, str], mode='t'):
        super().__init__(source=source, mode=mode, fmt='GTF')
    
    def parse(self, handle:IO[str]):
        """
        """
        records = self.iterate(handle)
        return records

    def iterate(self, handle:IO[str]):
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
                end=int(fields[4]),
                strand=strand
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

def parse(handle:Union[IO[str], str]):
    """
    """
    return GtfIterator(handle)

def read(handle:Union[IO[str], str]):
    """
    """
    iterator = parse(handle)
    try:
        record = next(iterator)
    except StopIteration:
        raise ValueError('No records found in handle')
    
    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    return record
    