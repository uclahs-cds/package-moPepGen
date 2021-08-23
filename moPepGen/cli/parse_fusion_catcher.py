""" Module for FusionCatcher parser """
from typing import List
from pathlib import Path
import argparse
from moPepGen import logger, seqvar, parser
from .common import add_args_reference, add_args_verbose, print_start_message,\
    print_help_if_missing_args, load_references, generate_metadata


# pylint: disable=W0212
def add_subparser_parse_fusion_catcher(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseFusionCatcher """

    p = subparsers.add_parser(
        name='parseFusionCatcher',
        help='Parse FusionCatcher result for moPepGen to call variant peptides.',
        description='Parse the FusionCatcher result to TVF format of variant'
        'records for moPepGen to call variant peptides. The genome'
    )
    p.add_argument(
        '-f', '--fusion',
        type=Path,
        help="Path to the FusionCatcher's output file.",
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
    add_args_reference(p)
    add_args_verbose(p)
    p.set_defaults(func=parse_fusion_catcher)
    print_help_if_missing_args(p)

def parse_fusion_catcher(args:argparse.Namespace) -> None:
    """ Parse FusionCatcher output and save it in TVF format. """
    # unpack args
    fusion = args.fusion
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    genome, anno, _ = load_references(args=args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    for record in parser.FusionCatcherParser.parse(fusion):
        var_records = record.convert_to_variant_records(anno, genome)
        variants.extend(var_records)

    if args.verbose:
        logger(f'FusionCatcher output {fusion} loaded.')

    variants.sort()

    if args.verbose:
        logger('Variants sorted.')

    metadata = generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    if args.verbose:
        logger("Variants written to disk.")
