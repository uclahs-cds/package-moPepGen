""" This module takes in a variant peptide FASTA file and encodes the fasta
headers with 36-digit long UUID (32 digits of hexdecimal characters and 4
hypens). The original headers together with  the UUIDs arer saved into the
dict file at the same location. This resolves the problems that some proteomic
search engines have strick requirement on the FASTA header length. """
from __future__ import annotations
import argparse
from pathlib import Path
from typing import TYPE_CHECKING, Dict
import uuid
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from moPepGen.cli import common


if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

INPUT_FILE_FORMATS = ['.fa', '.fasta']
OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_encode_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen encodeFasta """
    parser:argparse.ArgumentParser = subparser.add_parser(
        name='encodeFasta',
        help='Encode variant peptide FASTA file header.',
        description='Encode variant peptide FASTA file header.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=parser, formats=INPUT_FILE_FORMATS,
        message='Input FASTA file, must be generated by moPepGen callVariant.'
    )
    common.add_args_output_path(
        parser=parser, formats=OUTPUT_FILE_FORMATS
    )
    common.add_args_decoy(parser)
    common.add_args_debug_level(parser)
    common.print_help_if_missing_args(parser)
    parser.set_defaults(func=encode_fasta)
    return parser

def encode_fasta(args:argparse.Namespace) -> None:
    """ encode fasta """
    common.print_start_message(args)

    input_path:Path = args.input_path
    output_path:Path = args.output_path
    common.validate_file_format(
        input_path, INPUT_FILE_FORMATS, check_readable=True
    )
    common.validate_file_format(
        output_path, OUTPUT_FILE_FORMATS, check_writable=False
    )

    fasta_dict = output_path.with_suffix(output_path.suffix + '.dict')

    with open(input_path, 'rt') as in_handle, \
            open(output_path, 'wt') as out_handle, \
            open(fasta_dict, 'wt') as dict_handle:
        record2title = lambda x: x.description
        id_mapper:Dict[str,str] = {}
        writer = FastaIO.FastaWriter(out_handle, record2title=record2title)
        record:SeqRecord
        for record in SeqIO.parse(in_handle, 'fasta'):
            is_decoy = is_decoy_sequence(record.description, args.decoy_string,
                args.decoy_string_position)
            if is_decoy:
                header = get_real_header(record.description, args.decoy_string,
                    args.decoy_string_position)
            else:
                header = record.description

            if header in id_mapper:
                index = id_mapper[header]
            else:
                index = str(uuid.uuid4())
                dict_handle.write(f"{index}\t{header}" + '\n')
                id_mapper[header] = index

            if is_decoy:
                index = get_decoy_header(index, args.decoy_string,
                    args.decoy_string_position)

            record.description = index
            writer.write_record(record)

def is_decoy_sequence(header:str, decoy_string:str, decoy_string_position:str) -> bool:
    """ Checks if the sequence is decoy """
    if decoy_string_position == 'prefix':
        return header.startswith(decoy_string)
    if decoy_string_position == 'suffix':
        return header.endswith(decoy_string)
    raise ValueError(f"decoy_string_position of {decoy_string_position} is not supported.")

def get_real_header(header:str, decoy_string:str, decoy_string_position:str) -> str:
    """ Get the real header if it is a decoy sequence """
    if decoy_string_position == 'prefix':
        return header[len(decoy_string):]
    if decoy_string_position == 'suffix':
        return header[:-len(decoy_string)]
    raise ValueError(f"decoy_string_position of {decoy_string_position} is not supported.")

def get_decoy_header(header:str, decoy_string:str, decoy_string_position:str) -> str:
    """ Get the decoyed version of the header """
    if decoy_string_position == 'prefix':
        return decoy_string + header
    if decoy_string_position == 'suffix':
        return header + decoy_string
    raise ValueError(f"decoy_string_position of {decoy_string_position} is not supported.")
