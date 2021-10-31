""" Module for STAR-Fusion parser """
from __future__ import annotations
from typing import List, TYPE_CHECKING
from moPepGen import logger, seqvar, parser
from .common import add_args_reference, add_args_verbose, add_args_source,\
    print_start_message, print_help_if_missing_args, load_references, \
    generate_metadata


if TYPE_CHECKING:
    import argparse

# pylint: disable=W0212
def add_subparser_parse_star_fusion(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseSTARFusion """

    p = subparsers.add_parser(
        name='parseSTARFusion',
        help='Parse STAR-Fusion result for moPepGen to call variant peptides.',
        description='Parse the STAR-Fusion result to GVF format of variant'
        'records for moPepGen to call variant peptides.'
    )

    p.add_argument(
        '-f', '--fusion',
        type=str,
        help="Path to the STAR-Fusion's output file.",
        metavar='',
        required=True
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar='',
        required=True
    )
    p.add_argument(
        '--min-est-j',
        help='Minimal estimated junction reads. Defaults to 5.0',
        type=float,
        default=5.0,
        metavar=''
    )
    add_args_source(p)
    add_args_reference(p, proteome=False)
    add_args_verbose(p)
    p.set_defaults(func=parse_star_fusion)
    print_help_if_missing_args(p)
    return p

def parse_star_fusion(args:argparse.Namespace) -> None:
    """ Parse the STAR-Fusion's output and save it in GVF format. """
    # unpack args
    fusion = args.fusion
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    genome, anno, *_ = load_references(args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    for record in parser.STARFusionParser.parse(fusion):
        if record.est_j < args.min_est_j:
            continue
        var_records = record.convert_to_variant_records(anno, genome)
        variants.extend(var_records)

    if args.verbose:
        logger(f'STAR-Fusion output {fusion} loaded.')

    variants.sort()

    if args.verbose:
        logger('Variants sorted.')

    metadata = generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    if args.verbose:
        logger('Variant info written to disk.')
