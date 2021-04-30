from typing import IO, Union
from Bio.SeqIO.Interfaces import SequenceIterator
from moPepGen.GTFRecord import GTFRecord


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
            
            record = GTFRecord(
                seqname=fields[0],
                source=fields[1],
                feature=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                score=None if fields[5] == '.' else float(fields[5]),
                strand=fields[6],
                frame=None if fields[7] == '.' else float(fields[7]),
                attributes={key:val.strip('"') for key,val in \
                    [field.strip().split(' ') for field in \
                        fields[8].rstrip(';').split(';')]}
            )
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
    