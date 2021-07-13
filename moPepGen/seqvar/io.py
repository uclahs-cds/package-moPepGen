""" Module for moPepGen seqvar IO """
from typing import Iterable
import tempfile
from moPepGen.seqvar.TVFMetadata import TVFMetadata
from moPepGen.seqvar.VariantRecord import VariantRecord, ATTRS_START
from moPepGen.SeqFeature import FeatureLocation


def parse(path:str) -> Iterable[VariantRecord]:
    """ Parse a TVF file.

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
            start=int(fields[1]) - 1
            ref = fields[3]
            alt = fields[4]
            attrs = {}
            for field in fields[7].split(';'):
                key, val = field.split('=')
                val = val.strip('"')
                if key in ATTRS_START:
                    val = str(int(val) - 1)
                attrs[key] = val

            if not alt.startswith('<'):
                end = start + len(ref)
                _type = 'SNV' if len(ref) == len(alt) else 'INDEL'
            elif alt == '<FUSION>':
                end = start + 1
                _type = 'Fusion'
            elif alt == '<DEL>':
                end = int(attrs['END'])
                _type = 'Deletion'
            elif alt == '<INS>':
                end = start + 1
                _type = 'Insertion'
            elif alt == '<SUB>':
                end = int(attrs['END'])
                _type = 'Substitution'
            else:
                raise ValueError('Alt type not supported.')

            _id = fields[2]

            record = VariantRecord(
                location=FeatureLocation(
                    seqname=transcript_id,
                    start=start,
                    end=end
                ),
                ref=ref,
                alt=alt,
                _type=_type,
                _id=_id,
                attrs=attrs
            )
            yield record
            line = next(handle, None)


def write(variants:Iterable[VariantRecord], output_path:str,
        metadata:TVFMetadata):
    """ Write variants to an out handle.

    Args:
        variants (Iterable[seqvar.VariantRecord]): An iterable of variant
            records to write.
        handle (IO): The destination handle to write out.
        mode (str): If 'w', the header will be written, otherwise not.
    """
    temp_file = tempfile.TemporaryFile(mode='w+t')
    headers = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    temp_file.write('#' + '\t'.join(headers) + '\n')
    for record in variants:
        line = record.to_tvf()
        temp_file.write(line + '\n')
        metadata.add_info(record.type)

    out_file = open(output_path, 'w')

    for line in metadata.to_strings():
        out_file.write(line + '\n')

    temp_file.seek(0)
    for line in temp_file:
        out_file.write(line)

    out_file.close()
