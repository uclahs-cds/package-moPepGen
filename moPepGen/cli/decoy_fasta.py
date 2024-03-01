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
from moPepGen import aa, get_logger
from moPepGen.cli import common
from moPepGen.aa.expasy_rules import EXPASY_RULES


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
        '--enzyme',
        type=str,
        choices=[None, *EXPASY_RULES.keys()],
        help='Enzymatic cleavage rule. Amino acids at cleavage sites will be'
        ' kept unmodified. Set it to None to turn off this behavior.',
        metavar='<value>',
        default=None
    )
    parser.add_argument(
        '--non-shuffle-pattern',
        type=str,
        help='Residues to not shuffle and keep at the original position.'
        ' Separate by common (e.g. "K,R")',
        metavar='<value>',
        default=''
    )
    parser.add_argument(
        '--shuffle-max-attempts',
        type=int,
        help='Maximal attempts to shuffle a sequence to avoid any identical'
        ' decoy sequence.',
        metavar='<number>',
        default=30
    )
    parser.add_argument(
        '--keep-peptide-nterm',
        type=str,
        choices=['true', 'false'],
        metavar='<choice>',
        default='true',
        help='Whether to keep the peptide N terminus constant.'
    )
    parser.add_argument(
        '--keep-peptide-cterm',
        type=str,
        choices=['true', 'false'],
        default='true',
        metavar='<choice>',
        help='Whether to keep the peptide C terminus constant.'
    )
    parser.add_argument(
        '--seed',
        type=int,
        help='Random seed number.',
        metavar='<number>',
        default=None
    )
    parser.add_argument(
        '--order',
        type=str,
        choices=['juxtaposed', 'target_first', 'decoy_first'],
        help='Order of target and decoy sequences to write in the output FASTA.',
        metavar='<choice>',
        default='juxtaposed'
    )
    common.add_args_debug_level(parser)
    common.print_help_if_missing_args(parser)
    parser.set_defaults(func=decoy_fasta)
    return parser

class _Summary():
    """ Summary """
    def __init__(self, n_decoy:int=0, n_overlap:int=0):
        """ Constructor """
        self.n_decoy = n_decoy
        self.n_overlap = n_overlap

    def log_summary(self):
        """ Print summary to stdout """
        get_logger().info("Number of decoy sequences created: %s", self.n_decoy)
        get_logger().info(
            "Number of decoy sequences overlap with either target or decoy: %i",
            self.n_overlap
        )

class DecoyFasta():
    """ Decoy Fasta """
    def __init__(self, input_path:Path, output_path:Path, method:str,
            enzyme:str, keep_peptide_nterm:bool, keep_peptide_cterm:bool,
            non_shuffle_pattern:List[str], shuffle_max_attempts:int,  seed:int,
            decoy_string:str, decoy_string_position:str, order:str,
            target_db:List[SeqRecord]=None, _target_pool:Set[SeqRecord]=None,
            decoy_db:List[SeqRecord]=None, _decoy_pool:Set[SeqRecord]=None,
            _summary:_Summary=None):
        """ Constructor """
        self.input_path = input_path
        self.output_path = output_path
        self.method = method
        self.seed = seed
        self.enzyme = enzyme
        self.keep_peptide_nterm = keep_peptide_nterm
        self.keep_peptide_cterm = keep_peptide_cterm
        self.non_shuffle_pattern = non_shuffle_pattern
        self.shuffle_max_attempts = shuffle_max_attempts
        self.decoy_string = decoy_string
        self.decoy_string_position = decoy_string_position
        self.order = order
        self.target_db = target_db or []
        self._target_pool = _target_pool or set()
        self.decoy_db = decoy_db or []
        self._decoy_pool = _decoy_pool or set()
        self._summary = _summary or _Summary()

    @classmethod
    def from_args(cls, args:argparse.ArgumentParser):
        """ Create instance from args """
        common.print_start_message(args)

        input_path:Path = args.input_path
        output_path:Path = args.output_path
        common.validate_file_format(
            input_path, INPUT_FILE_FORMATS, check_readable=True
        )
        common.validate_file_format(
            output_path, OUTPUT_FILE_FORMATS, check_writable=True
        )

        keep_peptide_nterm = args.keep_peptide_nterm == 'true'
        keep_peptide_cterm = args.keep_peptide_cterm == 'true'

        non_shuffle_pattern = args.non_shuffle_pattern.split(',')

        return cls(
            input_path=input_path,
            output_path=output_path,
            method=args.method,
            enzyme=args.enzyme,
            keep_peptide_nterm=keep_peptide_nterm,
            keep_peptide_cterm=keep_peptide_cterm,
            non_shuffle_pattern=non_shuffle_pattern,
            shuffle_max_attempts=args.shuffle_max_attempts,
            seed=args.seed,
            decoy_string=args.decoy_string,
            decoy_string_position=args.decoy_string_position,
            order=args.order
        )

    def find_fixed_indices(self, seq:Seq) -> List[int]:
        """ Find the peptide indices to be fixed """
        fixed_indices = []

        if self.enzyme is not None:
            rule = self.enzyme
            exception = 'trypsin_expection' if self.enzyme == 'trypsin' else None
            fixed_indices += aa.AminoAcidSeqRecord(seq) \
                .find_all_enzymatic_cleave_sites(rule, exception)

        for i, it in enumerate(seq):
            if i == 0 and self.keep_peptide_nterm:
                fixed_indices.append(i)
                continue
            if i == len(seq) - 1 and self.keep_peptide_cterm:
                fixed_indices.append(i)
                continue
            if it in self.non_shuffle_pattern:
                fixed_indices.append(i)
        return fixed_indices

    @staticmethod
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

    @staticmethod
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

    def generate_decoy_sequence(self, seq:SeqRecord) -> SeqRecord:
        """ Generate decoy sequence """
        fixed_indices = self.find_fixed_indices(seq.seq)
        if self.method == 'reverse':
            decoy_seq = self.reverse_sequence(seq.seq, fixed_indices)
            if decoy_seq in self._target_pool or decoy_seq in self._decoy_pool:
                self._summary.n_overlap += 1
        elif self.method == 'shuffle':
            attempts = 0
            while True:
                decoy_seq = self.shuffle_sequence(seq.seq, fixed_indices)
                attempts += 1
                if decoy_seq not in self._target_pool and decoy_seq not in self._decoy_pool:
                    break
                if attempts >= self.shuffle_max_attempts:
                    self._summary.n_overlap += 1
                    break
        else:
            raise ValueError(f'Method {self.method} is not supported.')

        self._summary.n_decoy += 1

        self._decoy_pool.add(decoy_seq)

        if self.decoy_string_position == 'prefix':
            decoy_header = self.decoy_string + seq.description
        else:
            decoy_header = seq.description + self.decoy_string

        self.decoy_db.append(SeqRecord(decoy_seq, description=decoy_header))

    def iterate_target_decoy_database(self) -> Iterable[SeqRecord]:
        """ Iterate through target and decoy database with specified order """
        if self.order == 'juxtaposed':
            for i, target_seq in enumerate(self.target_db):
                yield target_seq
                yield self.decoy_db[i]
        elif self.order == 'target_first':
            for seq in self.target_db:
                yield seq
            for seq in self.decoy_db:
                yield seq
        elif self.order == 'decoy_first':
            for seq in self.decoy_db:
                yield seq
            for seq in self.target_db:
                yield seq
        else:
            raise ValueError(f'Order {self.order} is not supported.')

    def write(self):
        """ Write decoy database to fasta """
        with open(self.output_path, 'wt') as handle:
            record2title = lambda x: x.description
            writer = FastaIO.FastaWriter(handle, record2title=record2title)
            for seq in self.iterate_target_decoy_database():
                writer.write_record(seq)

    def main(self):
        """ Create decoy database """
        logger = get_logger()
        with open(self.input_path, 'rt') as handle:
            self.target_db = list(SeqIO.parse(handle, format='fasta'))
            self.target_db.sort(key=lambda x: x.seq)

        self._target_pool = {x.seq for x in self.target_db}

        logger.info('Input database FASTA file loaded.')

        if self.seed is not None:
            random.seed(self.seed)

        for seq in self.target_db:
            self.generate_decoy_sequence(seq)

        logger.info('Decoy sequences created.')

        self._summary.log_summary()

        self.write()
        logger.info('Decoy database written.')

def decoy_fasta(args:argparse.ArgumentParser):
    """ Generate decoy database fasta file for library searching. """
    DecoyFasta.from_args(args).main()
