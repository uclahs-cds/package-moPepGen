""" Module for moPepGen seqvar IO """
from typing import Iterable, IO
from moPepGen import seqvar
from moPepGen.SeqFeature import FeatureLocation


def parse(path:str) -> Iterable[seqvar.VariantRecord]:
    """ Parse a moPepGen variant file.

    Args:
        path (str): Path to the moPepGen variant file.
        type (str):
    Return:
        An iterable of seq
    """
    with open(path, 'r') as handle:
        line = next(handle, None)
        while line:
            if line.startswith('#'):
                line = next(handle, None)
                continue
            fields = line.rstrip().split('\t')
            transcript_id = fields[0]
            record = seqvar.VariantRecord(
                location=FeatureLocation(
                    seqname=transcript_id,
                    start=int(fields[1]),
                    end=int(fields[2])
                ),
                ref=fields[3],
                alt=fields[4],
                _type=fields[5],
                _id=fields[6]
            )
            yield record
            line = next(handle, None)


def write(variants:Iterable[seqvar.VariantRecord], handle:IO, mode:str='w'):
    """"""
    if mode == 'w':
        headers = ['transcript_id', 'start', 'end', 'ref', 'alt', 'type', 'id']
        handle.write('#' + '\t'.join(headers) + '\n')
    for record in variants:
        record:seqvar.VariantRecord
        transcript_id = record.location.seqname
        line = [transcript_id, str(int(record.location.start)),
            str(int(record.location.end)), str(record.ref),
            str(record.alt), record.type, record.id]
        handle.write('\t'.join(line) + '\n')

