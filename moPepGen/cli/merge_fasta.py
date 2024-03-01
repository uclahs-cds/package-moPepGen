""" `MergeFasta` merges a serious of database FASTA files into one. This is
useful when working with multiplexed proteomic experiments such as TMT """
from __future__ import annotations
import argparse
from moPepGen import get_logger
from moPepGen.aa.VariantPeptidePool import VariantPeptidePool
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.fasta', '.fa']
OUTPUT_FILE_FORMATS = ['.fasta', '.fa']

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
    output_file = args.output_path
    common.validate_file_format(
        output_file, OUTPUT_FILE_FORMATS, check_writable=True
    )

    pool = None

    for file in input_files:
        with open(file, 'rt') as handle:
            if pool is None:
                pool = VariantPeptidePool.load(handle)
            else:
                second_pool = VariantPeptidePool.load(handle)
                for peptide in second_pool.peptides:
                    pool.add_peptide(
                        peptide=peptide,
                        canonical_peptides=set(),
                        skip_checking=True
                    )
            logger.info("Database FASTA file loaded: %s", file)

    if args.dedup_header:
        pool.remove_redundant_headers()
    pool.write(output_file)

    logger.info("Merged FASTA file saved to %s", output_file)
