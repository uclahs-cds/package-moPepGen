""" `indexGVF`  """
from __future__ import annotations
import argparse
from pathlib import Path
import tempfile
from moPepGen import logger, check_sha512
from moPepGen.seqvar import GVFIndex
from moPepGen.seqvar import GVFMetadata
from .common import add_args_quiet, print_help_if_missing_args, print_start_message


# pylint: disable=W0212
def add_subparser_index_gvf(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen indexGVF """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='indexGVF',
        help='Generate index for GVF file.',
        description='Generate index for GVF file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-i', '--input-gvf',
        type=Path,
        help='Input variant GVF file.',
        metavar='<file>',
        required=True
    )
    add_args_quiet(p)
    p.set_defaults(func=index_gvf)
    print_help_if_missing_args(p)
    return p


def index_gvf(args:argparse.Namespace):
    """ Generate GVF index """
    input_file:Path = args.input_gvf
    output_file = input_file.with_suffix(input_file.suffix + '.idx')
    with open(input_file, 'rb') as handle:
        sum_val = check_sha512(handle)

    print_start_message(args)

    with open(input_file, 'rt') as in_handle:
        metadata = GVFMetadata.parse(in_handle)
        is_circ_rna = metadata.is_circ_rna()
        with tempfile.TemporaryFile(mode='w+t') as temp_file:
            temp_file.write(f"# CHECKSUM={sum_val}\n")
            in_handle.seek(0)
            it = GVFIndex.iterate_pointer(handle=in_handle, is_circ_rna=is_circ_rna)
            for pointer in it:
                temp_file.write(pointer.to_line() + '\n')
            temp_file.seek(0)
            with open(output_file, 'wt') as out_handle:
                for line in temp_file:
                    out_handle.write(line)

    logger('GVF index generated.')
