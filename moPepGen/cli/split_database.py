""" `splitDatabase` takes the FASTA file with variant peptide sequences called
by [`callVariant`](./call-variant.md) with or without noncoding novel peptides
called by [`callNoncoding`](./call-noncoding.md), and splits peptide sequences
into databases. The split database FASTA files can be used for sequential
library searching. """
from __future__ import annotations
import argparse
from pathlib import Path
from moPepGen.aa import PeptidePoolSplitter
from moPepGen import SPLIT_DATABASE_KEY_SEPARATER, logger
from .common import add_args_reference, add_args_quiet, print_start_message,\
    print_help_if_missing_args, load_references


# pylint: disable=W0212
def add_subparser_split_database(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen splitDatabase """
    p = subparser.add_parser(
        name='splitDatabase',
        help='Split variant peptide FASTA database generated by moPepGen.',
        description='Split variant peptide FASTA database generated by'
        'moPepGen into separate files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '-i', '--input-variant',
        type=Path,
        help='Path to the input variant GVF files',
        metavar='<file>',
        nargs='+'
    )
    p.add_argument(
        '-r', '--variant-peptides',
        type=Path,
        help='Path to the variant peptide FASTA database file.',
        metavar='<file>'
    )
    p.add_argument(
        '-n', '--noncoding-peptides',
        type=Path,
        help='Pth the the noncoding peptide FASTA databse file.',
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--order-source',
        type=str,
        help='Order of sources, separate by comma. E.g., SNP,SNV,Fusion',
        metavar='<value>'
    )
    p.add_argument(
        '--group-source',
        type=str,
        help='Group sources. E.g., PointMutation:gSNP,sSNV INDEL:gINDEL,sINDEL',
        metavar='<value>',
        nargs='*'
    )
    p.add_argument(
        '--max-source-groups',
        type=int,
        help='Maximal number of different source groups to be separate into'
        'individual database FASTA files. Defaults to 1',
        default=1
    )
    p.add_argument(
        '--additional-split',
        type=str,
        help='For peptides that were not already split into FASTAs up to'
        'max_source_groups, those involving the following source will be split'
        'into additional FASTAs with decreasing priority',
        metavar='<value>',
        default=None
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=Path,
        help='Output prefix',
        metavar = '<value>'
    )

    add_args_reference(p, genome=False, proteome=False)
    add_args_quiet(p)
    print_help_if_missing_args(p)
    p.set_defaults(func=split_database)
    return p

def split_database(args:argparse.Namespace) -> None:
    """ Split peptide database """
    print_start_message(args)

    _, anno, *_ = load_references(args, load_genome=False, \
        load_proteome=False, load_canonical_peptides=False)

    source_order = {val:i for i,val in  enumerate(args.order_source.split(','))}\
        if args.order_source else None

    group_map = None
    if args.group_source:
        group_map = {}
        for it in args.group_source:
            key, val = it.split(':')
            group_map[key] = val.split(',')

    splitter = PeptidePoolSplitter(order=source_order, group_map=group_map)

    with open(args.variant_peptides, 'rt') as handle:
        splitter.load_database(handle)

    if args.noncoding_peptides:
        with open(args.noncoding_peptides) as handle:
            splitter.load_database(handle)

    for file in args.variant_gvf:
        with open(file, 'rt') as handle:
            splitter.load_gvf(handle)

    additional_split = args.additional_split or []
    sep = SPLIT_DATABASE_KEY_SEPARATER
    additional_split = [{x.split(sep)} for x in additional_split]
    splitter.split(args.max_source_groups, additional_split, anno)

    if not args.quiet:
        logger('Database split finished')

    splitter.write(args.output_prefix)

    if not args.quiet:
        logger('Split databases saved to disk.')
