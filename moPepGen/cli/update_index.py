""" `updateIndex` creates new canoincal peptide pool to a moPepGen index dir
created by `generateIndex`. This command is useful when you want to use a
different enzyme or cleavage settings (*e.g.*, `miscleavages`, `min_length`) with
the same reference genome, annotation, and proteome. """
from __future__ import annotations
import argparse
import sys
from pathlib import Path
from moPepGen import params, get_logger
from moPepGen.index import IndexDir
from moPepGen.cli import common


# pylint: disable=W0212
def add_subparser_update_index(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen generateIndex """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='updateIndex',
        help='Update moPepGen index',
        description='Update moPepGen index',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '--index-dir',
        type=Path,
        help='Path to the directory of index files generated by moPepGen'
        ' generateIndex. If given, --genome-fasta, --proteome-fasta and'
        ' --anntotation-gtf will be ignored.',
        metavar='<file>',
        nargs='?',
        default=None
    )
    p.add_argument(
        '-f', '--force',
        action='store_true',
        help='Force write data to index dir.'
    )
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)
    p.set_defaults(func=update_index)
    common.print_help_if_missing_args(p)
    return p

def update_index(args:argparse.Namespace):
    """ update index """
    logger = get_logger()

    common.print_start_message(args)

    index_dir = IndexDir(args.index_dir)
    index_dir.validate_metadata()

    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    min_length:int = int(args.min_length)
    max_length:int = int(args.max_length)
    exception = args.cleavage_exception

    cleavage_params = params.CleavageParams(
        enzyme=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    pool_exists = index_dir.metadata.get_canonical_pool(cleavage_params) is not None
    if pool_exists and not args.force:
        logger.error(
            'Canonical peptide pool already exists with the given parameters.'
        )
        sys.exit(1)

    # load reference data
    anno = index_dir.load_annotation()
    logger.info('Genomic annotation loaded.')
    proteome = index_dir.load_proteome()
    logger.info('Proteome loaded.')

    # create canonical peptide pool
    canonical_peptides = proteome.create_unique_peptide_pool(
        anno=anno, rule=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    logger.info('Canoincal peptide pool generated.')
    index_dir.save_canonical_peptides(canonical_peptides, cleavage_params, override=args.force)
    logger.info('Canoincal peptide pool saved.')

    if not pool_exists:
        index_dir.save_metadata()
        logger.info('metadata.json updated.')
