""" Module for rMATS parser """
from __future__ import annotations
from typing import Dict, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen import logger, seqvar
from moPepGen.parser import RMATSParser
from .common import add_args_reference, add_args_verbose, print_start_message,\
    print_help_if_missing_args, load_references, generate_metadata


if TYPE_CHECKING:
    import argparse

# pylint: disable=W0212
def add_subparser_parse_rmats(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseRMATs """
    ## parser_rmats
    p = subparsers.add_parser(
        name='parseRMATS',
        help='Parse rMATS result for moPepGen to call variant peptides.',
        description='Parse the rMATS result to GVF format of variant'
        'records for moPepGen to call variant peptides.'
    )

    p.add_argument(
        '--skipped-exon',
        type=Path,
        help="Skipped exon junction count txt file.",
        metavar='',
        default=None,
        dest='skipped_exon'
    )
    p.add_argument(
        '--alternative-5-splicing',
        type=Path,
        help="Alternative 5' splicing junction count txt file.",
        metavar='',
        default=None,
        dest='alternative_5_splicing'
    )
    p.add_argument(
        '--alternative-3-splicing',
        type=Path,
        help="Alternative 3' splicing junction count txt file.",
        metavar='',
        default=None,
        dest='alternative_3_splicing'
    )
    p.add_argument(
        '--mutually-exclusive-exons',
        type=Path,
        help="Mutually exclusive junction count txt file.",
        metavar='',
        default=None,
        dest='mutually_exclusive_exons'
    )
    p.add_argument(
        '--retained-intron',
        type=Path,
        help="Retained intron junction count txt file.",
        metavar='',
        default=None,
        dest='retained_intron'
    )

    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar=''
    )

    add_args_reference(p, proteome=False)
    add_args_verbose(p)
    p.set_defaults(func=parse_rmats)
    print_help_if_missing_args(p)

def parse_rmats(args:argparse.Namespace) -> None:
    """ Parse rMATS results into TSV """
    skipped_exon = args.skipped_exon
    alternative_5 = args.alternative_5_splicing
    alternative_3 = args.alternative_3_splicing
    mutually_exclusive = args.mutually_exclusive_exons
    retained_intron = args.retained_intron
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.tvf'

    print_start_message(args)

    genome, anno, *_ = load_references(args, load_canonical_peptides=False)

    variants:Dict[str,Set[seqvar.VariantRecord]] = {}
    rmats_outputs = [
        ('SE', skipped_exon), ('A5SS', alternative_5), ('A3SS', alternative_3),
        ('MXE', mutually_exclusive), ('RI', retained_intron)
    ]
    for event_type, path in rmats_outputs:
        if path:
            for record in RMATSParser.parse(path, event_type):
                var_records = record.convert_to_variant_records(anno, genome)
                for var_record in var_records:
                    transcript_id = var_record.location.seqname
                    if transcript_id not in variants:
                        variants[transcript_id] = set()
                    variants[transcript_id].add(var_record)

    variants_sorted = []
    for val in variants.values():
        val = list(val)
        val.sort()
        variants_sorted.extend(val)

    if args.verbose:
        logger('Variants sorted.')

    metadata = generate_metadata(args)
    seqvar.io.write(variants_sorted, output_path, metadata)

    if args.verbose:
        logger("Variants written to disk.")


if __name__ == '__main__':
    data = argparse.Namespace()
    data.skipped_exon = Path('test/files/alternative_splicing/rmats_se.txt')
    data.alternative_5_splicing = Path('test/files/alternative_splicing/rmats_a5ss.txt')
    data.alternative_3_splicing = Path('test/files/alternative_splicing/rmats_a3ss.txt')
    data.mutually_exclusive_exons = Path('test/files/alternative_splicing/rmats_mxe.txt')
    data.retained_intron = Path('test/files/alternative_splicing/rmats_ri.txt')
    data.index_dir = Path('test/files/index')
    data.output_prefix = 'test/files/alternative_splicing/rmats'
    data.verbose = True
    parse_rmats(data)
