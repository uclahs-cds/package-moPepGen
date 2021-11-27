""" `parseSTARFusion` takes the identified fusion transcript results from
[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) and save as a
GVF file. The GVF file can be later used to call variant peptides using
[callVariant](call-variant.md)."""
from __future__ import annotations
import argparse
from typing import List
from moPepGen import logger, seqvar, parser
from .common import add_args_output_prefix, add_args_reference, \
    add_args_quiet, add_args_source, print_start_message, \
    print_help_if_missing_args, load_references, generate_metadata


# pylint: disable=W0212
def add_subparser_parse_star_fusion(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseSTARFusion """

    p = subparsers.add_parser(
        name='parseSTARFusion',
        help='Parse STAR-Fusion result for moPepGen to call variant peptides.',
        description='Parse STAR-Fusion output to GVF format of variant'
        ' records for moPepGen to call variant peptides.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '-f', '--fusion',
        type=str,
        help="Path to STAR-Fusion's output file.",
        metavar='<file>',
        required=True
    )
    add_args_output_prefix(p)
    p.add_argument(
        '--min-est-j',
        help='Minimal estimated junction reads to be included.',
        type=float,
        default=5.0,
        metavar='<number>'
    )
    add_args_source(p)
    add_args_reference(p, proteome=False)
    add_args_quiet(p)
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

    if not args.quiet:
        logger(f'STAR-Fusion output {fusion} loaded.')

    variants.sort()

    if not args.quiet:
        logger('Variants sorted.')

    metadata = generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    if not args.quiet:
        logger('Variant info written to disk.')
