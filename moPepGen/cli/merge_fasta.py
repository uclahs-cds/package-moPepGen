""" `MergeFasta` merges a serious of database FASTA files into one. This is
useful when working with multiplexed proteomic experiments such as TMT """
from __future__ import annotations
from typing import TYPE_CHECKING
import argparse
import os
from pathlib import Path
import pickle
from Bio import SeqIO
from moPepGen import get_logger
from moPepGen.aa.VariantPeptidePool import VariantPeptidePool
from moPepGen.svgraph.VariantPeptideTable import (
    VariantPeptideTable,
    get_peptide_table_path,
    get_peptide_table_path_temp
)
from moPepGen.cli import common


if TYPE_CHECKING:
    from typing import Iterable, Set
    from Bio.Seq import Seq
    from logging import Logger

INPUT_FILE_FORMATS = ['.fasta', '.fa']
OUTPUT_FILE_FORMATS = ['.fasta', '.fa']
DENYLIST_FORMATS = ['.fasta', '.fa', '.pkl']

# pylint: disable=W0212
def add_subparser_merge_fasta(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen MergeFasta """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='mergeFasta',
        help='Merge multiple variant peptide FASTA files into one.',
        description='Merge multiple variant peptide FASTA files into one.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=parser, formats=INPUT_FILE_FORMATS, plural=True,
        message='File path to the variant peptide FASTA files.'
    )
    common.add_args_output_path(
        parser=parser, formats=OUTPUT_FILE_FORMATS
    )
    parser.add_argument(
        '--dedup-header',
        action='store_true',
        help='Remove duplicate FASTA header entries after merging.'
    )
    parser.add_argument(
        '--denylist',
        type=Path,
        nargs='*',
        default=[],
        help=f"Denylist of peptide sequences. Valid formats: {DENYLIST_FORMATS}"
    )
    common.add_args_debug_level(parser)
    parser.set_defaults(func=merge_fasta)
    common.print_help_if_missing_args(parser)
    return parser

def merge_fasta(args:argparse.Namespace):
    """ Merge mulitple variant peptide FASTA files into one. """
    logger = get_logger()
    input_files = args.input_path
    for file in input_files:
        common.validate_file_format(
            file, INPUT_FILE_FORMATS, check_readable=True
        )
    output_file:Path = args.output_path
    common.validate_file_format(
        output_file, OUTPUT_FILE_FORMATS, check_writable=True
    )

    for path in args.denylist:
        common.validate_file_format(
            path, DENYLIST_FORMATS, check_readable=True
        )

    denylist:Set[Seq] = set()
    path:Path
    for path in args.denylist:
        if path.suffix in ['.fasta', 'fa']:
            for rec in SeqIO.parse(path, 'fasta'):
                denylist.add(rec.seq)
        elif path.suffix in ['.pkl']:
            with open(path, 'rb') as handle:
                peptides = pickle.load(handle)
                for peptide in peptides:
                    denylist.add(peptide)

    if all_fasta_have_table(input_files):
        temp_file = get_peptide_table_path_temp(output_file)
        table_file = get_peptide_table_path(output_file)
        with open(temp_file, 'w+') as handle:
            peptide_table = VariantPeptideTable(handle=handle)
            peptide_table.write_header()
            merge_peptide_table(
                files=input_files,
                peptide_table=peptide_table,
                denylist=denylist,
                logger=logger
            )
            peptide_table.write_fasta(output_file)
            peptide_table.sort_table(table_file)
        os.remove(temp_file)
    else:
        pool = merge_peptide_fasta(
            files=input_files,
            denylist=denylist,
            logger=logger
        )
        if args.dedup_header:
            pool.remove_redundant_headers()
        pool.write(output_file)

    logger.info("Merged FASTA file saved to %s", output_file)

def all_fasta_have_table(files:Iterable[Path]):
    """ Check wether all fasta files have the peptide table """
    return all(get_peptide_table_path(Path(str(path))).exists() for path in files)

def merge_peptide_fasta(files:Iterable[Path], denylist:Set[Seq], logger:Logger=None):
    """ Merge peptides from FASTA files. """
    pool = None

    for path in files:
        with open(path, 'rt') as handle:
            if pool is None:
                pool = VariantPeptidePool.load(handle)
                if denylist:
                    pool = pool.filter_denylist(denylist)
            else:
                second_pool = VariantPeptidePool.load(handle)
                if denylist:
                    second_pool = second_pool.filter_denylist(denylist)
                for peptide in second_pool.peptides:
                    pool.add_peptide(
                        peptide=peptide,
                        canonical_peptides=set(),
                        skip_checking=True
                    )
            if logger:
                logger.info("Database FASTA file loaded: %s", path)

    return pool

def merge_peptide_table(files:Iterable[Path], peptide_table:VariantPeptideTable,
        denylist:Set[Seq], logger:Logger=None):
    """ Merge peptides from FASTA table """
    for path in files:
        table_path = get_peptide_table_path(path)
        with open(table_path, 'rt') as handle:
            table = VariantPeptideTable(handle)
            table.generate_index()
            for peptide in table.index:
                if peptide in denylist:
                    continue
                annos = table.load_peptide_annotation(peptide)
                for anno in annos.values():
                    peptide_table.add_peptide(seq=peptide, peptide_anno=anno)
        if logger:
            logger.info("Database FASTA file loaded: %s", path)
    return peptide_table
