""" This module takes a FASTA file and creates a decoy database by shuffling or
reversing each sequence. The generated decoy database FASTA file can then be
used for library searching with proteomics data. """
from __future__ import annotations
import argparse
from pathlib import Path
import random
from typing import Iterable, List, Set
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from moPepGen import logger
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.fa', '.fasta']
OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_decoy_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen decoyFasta """
    parser:argparse.ArgumentParser = subparser.add_parser(
        name='decoyFasta',
        help='Generate decoy database FASTA file.',
        description='Generate decoy database FASTA file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=parser, formats=INPUT_FILE_FORMATS,
        message='Input FASTA file.'
    )
    common.add_args_output_path(
        parser=parser, formats=OUTPUT_FILE_FORMATS
    )
    common.add_args_decoy(parser)
    parser.add_argument(
        '--method',
        type=str,
        choices=['reverse', 'shuffle'],
        help='Method to be used to generate the decoy sequences from target'
        ' sequences.',
        default='reverse'
    )
    parser.add_argument(
        '--non-shuffle-pattern',
        type=str,
        help='Residues to not shuffle and keep at the original position.'
        ' Separate by common (e.g. "K,R")',
        default=''
    )
    parser.add_argument(
        '--keep-peptide-nterm',
        type=str,
        choices=['true', 'false'],
        default='true',
        help='Whether to keep the peptide N terminus constant.'
    )
    parser.add_argument(
        '--keep-peptide-cterm',
        type=str,
        choices=['true', 'false'],
        default='true',
        help='Whether to keep the peptide C terminus constant.'
    )
    parser.add_argument(
        '--seed',
        type=int,
        help='Random seed number.',
        default=None
    )
    parser.add_argument(
        '--order',
        type=str,
        choices=['juxtaposed', 'target_first', 'decoy_first'],
        help='Order of target and decoy sequences to write in the output FASTA.',
        default='juxtaposed'
    )
    common.add_args_quiet(parser)
    common.print_help_if_missing_args(parser)
    parser.set_defaults(func=decoy_fasta)


def decoy_fasta(args:argparse.Namespace):
    """ Create decoy database """
    common.print_start_message(args)

    input_path:Path = args.input_path
    output_path:Path = args.output_path
    common.validate_file_format(input_path, INPUT_FILE_FORMATS, True)
    common.validate_file_format(input_path, OUTPUT_FILE_FORMATS, True)

    keep_peptide_nterm = args.keep_peptide_nterm == 'true'
    keep_peptide_cterm = args.keep_peptide_cterm == 'true'

    non_shuffle_pattern = args.non_shuffle_pattern.split(',')

    with open(input_path, 'rt') as handle:
        target_seqs:List[SeqRecord] = list(SeqIO.parse(handle, format='fasta'))

    target_pool = {x.seq for x in target_seqs}

    if not args.quiet:
        logger('Input database FASTA file loaded.')

    if args.seed is not None:
        random.seed(args.seed)

    decoy_seqs = []
    decoy_pool = set()
    for seq in target_seqs:
        decoy_seq = generate_decoy_sequence(
            seq, method=args.method, decoy_string=args.decoy_string,
            decoy_string_position=args.decoy_string_position,
            keep_nterm=keep_peptide_nterm,
            keep_cterm=keep_peptide_cterm,
            non_shuffle_pattern=non_shuffle_pattern,
            target_db=target_pool, decoy_db=decoy_pool
        )
        decoy_seqs.append(decoy_seq)
        decoy_pool.add(decoy_seq.seq)

    if not args.quiet:
        logger('Decoy sequences created.')

    with open(output_path, 'wt') as handle:
        record2title = lambda x: x.description
        writer = FastaIO.FastaWriter(handle, record2title=record2title)
        for seq in iterate_target_decoy_database(target_seqs, decoy_seqs, args.order):
            writer.write_record(seq)

    if not args.quiet:
        logger('Decoy database written.')

def generate_decoy_sequence(seq:SeqRecord, method:str, decoy_string:str,
        decoy_string_position:str, keep_nterm:bool, keep_cterm:bool,
        non_shuffle_pattern:List[str], target_db:Set[Seq], decoy_db:Set[Seq]
        ) -> SeqRecord:
    """ Generate decoy sequence """
    fixed_indices = find_fixed_indices(seq.seq, keep_nterm, keep_cterm, non_shuffle_pattern)
    if method == 'reverse':
        decoy_seq = reverse_sequence(seq.seq, fixed_indices)
    elif method == 'shuffle':
        while True:
            decoy_seq = shuffle_sequence(seq.seq, fixed_indices)
            if decoy_seq not in target_db and decoy_seq not in decoy_db:
                break
    else:
        raise ValueError(f'Method {method} is not supported.')

    if decoy_string_position == 'prefix':
        decoy_header = decoy_string + seq.description
    else:
        decoy_header = seq.description + decoy_string

    return SeqRecord(decoy_seq, description=decoy_header)

def reverse_sequence(seq:Seq, fixed_indices) -> Seq:
    """ Create a reversed sequence """
    seq = str(seq)
    reversed_indices = list(reversed([i for i,_ in enumerate(seq) if i not in fixed_indices]))
    shuffled_seq = []
    offset = 0
    i = 0
    while i < len(reversed_indices):
        if i + offset in fixed_indices:
            shuffled_seq.append(seq[i + offset])
            offset += 1
            continue
        shuffled_seq.append(seq[reversed_indices[i]])
        i += 1

    if len(shuffled_seq) < len(seq):
        shuffled_seq += list(seq[len(shuffled_seq) - len(seq):])
    return Seq(''.join(shuffled_seq))

def shuffle_sequence(seq:Seq, fixed_indices:List[int]) -> Seq:
    """ Create a shuffled sequence """
    seq = str(seq)
    indices_to_shuffle = [i for i,_ in enumerate(seq) if i not in fixed_indices]
    shuffled_indices = random.sample(indices_to_shuffle, len(indices_to_shuffle))
    shuffled_seq = []
    offset = 0
    i = 0
    while i < len(shuffled_indices):
        if i + offset in fixed_indices:
            shuffled_seq.append(seq[i + offset])
            offset += 1
            continue
        shuffled_seq.append(seq[shuffled_indices[i]])
        i += 1

    if len(shuffled_seq) < len(seq):
        shuffled_seq += list(seq[len(shuffled_seq) - len(seq):])
    return Seq(''.join(shuffled_seq))

def find_fixed_indices(seq:Seq, keep_nterm:bool, keep_cterm:bool, non_shuffle_pattern:List[str]
        ) -> List[int]:
    """ Find the peptide indices to be fixed """
    fixed_indices = []
    for i, it in enumerate(seq):
        if i == 0 and keep_nterm:
            fixed_indices.append(i)
            continue
        if i == len(seq) - 1 and keep_cterm:
            fixed_indices.append(i)
            continue
        if it in non_shuffle_pattern:
            fixed_indices.append(i)
    return fixed_indices

def iterate_target_decoy_database(target:List[SeqRecord], decoy:List[SeqRecord],
        order=str) -> Iterable[SeqRecord]:
    """ Iterate through target and decoy database with specified order """
    if order == 'juxtaposed':
        for i, target_seq in enumerate(target):
            yield target_seq
            yield decoy[i]
    elif order == 'target_first':
        for seq in target:
            yield seq
        for seq in decoy:
            yield seq
    elif order == 'decoy_first':
        for seq in decoy:
            yield seq
        for seq in target:
            yield seq
    else:
        raise ValueError(f'Order {order} isnot supported.')
