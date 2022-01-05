""" Module for moPepGen seqvar IO """
from typing import Iterable, Union, IO
from pathlib import Path
import tempfile
from moPepGen import GVF_HEADER
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen.seqvar.VariantRecord import VariantRecord, ATTRS_POSITION
from moPepGen.SeqFeature import FeatureLocation


T = Union[str, IO, Path]

def parse(handle:T) -> Iterable[VariantRecord]:
    """ Parse a GVF file.

    Args:
        path (str): Path to the moPepGen variant file.
        type (str):
    Return:
        An iterable of seq
    """
    if isinstance(handle, (str, Path)):
        with open(handle, 'r') as stream:
            for record in iterate(stream):
                yield record
    else:
        for record in iterate(handle):
            yield record

def iterate(handle:IO) -> Iterable[VariantRecord]:
    """ Iterate through a file-handle of a GVF file and returns an iterator
    of VariantRecord objects.

    Args:
        handle (IO): File-handle like object. Usually the TextIOWrapper
            returned by the open statement.
        it (str): An iterator if the iterating already starts.

    Example:
        >>> with open('test.gvf', 'rt') as handle:
                for record in iterate(handle):
                    print(record)
                    break
        <VariantRecord>
    """
    for line in handle:
        if line.startswith('#'):
            continue
        yield line_to_variant_record(line)


def line_to_variant_record(line:str) -> VariantRecord:
    """ """
    fields = line.rstrip().split('\t')
    gene_id = fields[0]
    start=int(fields[1]) - 1
    ref = fields[3]
    alt = fields[4]
    attrs = {}
    for field in fields[7].split(';'):
        key, val = field.split('=')
        val = val.strip('"')
        if key in ATTRS_POSITION:
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

    return VariantRecord(
        location=FeatureLocation(
            seqname=gene_id,
            start=start,
            end=end
        ),
        ref=ref,
        alt=alt,
        _type=_type,
        _id=_id,
        attrs=attrs
    )

def write(variants:Iterable[VariantRecord], output_path:str,
        metadata:GVFMetadata):
    """ Write variants to an out handle.

    Args:
        variants (Iterable[seqvar.VariantRecord]): An iterable of variant
            records to write.
        handle (IO): The destination handle to write out.
        mode (str): If 'w', the header will be written, otherwise not.
    """
    with tempfile.TemporaryFile(mode='w+t') as temp_file:
        temp_file.write('#' + '\t'.join(GVF_HEADER) + '\n')
        for record in variants:
            line = record.to_string()
            temp_file.write(line + '\n')
            metadata.add_info(record.type)
        with open(output_path, 'w') as out_file:
            for line in metadata.to_strings():
                out_file.write(line + '\n')

            temp_file.seek(0)
            for line in temp_file:
                out_file.write(line)
