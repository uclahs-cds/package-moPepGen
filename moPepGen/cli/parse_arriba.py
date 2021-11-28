""" `parseArriba` takes the identified fusion transcript results from
[Arriba](https://github.com/suhrig/arriba) and saves as a GVF file. The GVF
file can be later used to call variant peptides using
[callVariant](call-variant.md)."""
from typing import List
from pathlib import Path
import argparse
from moPepGen import logger, seqvar, parser, err
from .common import add_args_reference, add_args_verbose, add_args_source,\
    add_args_output_prefix, print_start_message,print_help_if_missing_args,\
    load_references, generate_metadata


# pylint: disable=W0212
def add_subparser_parse_arriba(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseArriba """

    p = subparsers.add_parser(
        name='parseArriba',
        help='Parse Arriba result for moPepGen to call variant peptides.',
        description='Parse the Arriba result to GVF format of variant'
        'records for moPepGen to call variant peptides.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-f', '--fusion',
        type=Path,
        help="Path to Arriba's output file.",
        metavar='<file>',
        required=True
    )
    add_args_output_prefix(p)
    add_args_source(p)
    add_args_reference(p, proteome=False)
    add_args_verbose(p)
    p.set_defaults(func=parse_arriba)
    print_help_if_missing_args(p)
    return p

def parse_arriba(args:argparse.Namespace) -> None:
    """ Parse Arriba output and save it in GVF format. """
    # unpack args
    fusion = args.fusion
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    genome, anno, *_ = load_references(args=args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    with open(fusion, 'rt') as handle:
        for record in parser.ArribaParser.parse(handle):
            if record.transcript_on_antisense_strand(anno):
                continue
            try:
                var_records = record.convert_to_variant_records(anno, genome)
            except err.GeneNotFoundError:
                continue
            variants.extend(var_records)

    if args.verbose:
        logger(f'Arriba output {fusion} loaded.')

    variants.sort()

    if args.verbose:
        logger('Variants sorted.')

    metadata = generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    if args.verbose:
        logger("Variants written to disk.")