""" `parseSTARFusion` takes the identified fusion transcript results from
[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) and save as a
GVF file. The GVF file can be later used to call variant peptides using
[callVariant](call-variant.md)."""
from __future__ import annotations
import argparse
from typing import List
from moPepGen import get_logger, seqvar, parser, err
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.tsv', '.txt']
OUTPUT_FILE_FORMATS = ['.gvf']

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

    common.add_args_input_path(
        parser=p,
        formats=INPUT_FILE_FORMATS,
        message="File path to STAR-Fusion's output file."
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--min-est-j',
        help='Minimal estimated junction reads to be included.',
        type=float,
        default=5.0,
        metavar='<number>'
    )
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_star_fusion)
    common.print_help_if_missing_args(p)
    return p

def parse_star_fusion(args:argparse.Namespace) -> None:
    """ Parse the STAR-Fusion's output and save it in GVF format. """
    logger = get_logger()
    # unpack args
    fusion = args.input_path
    common.validate_file_format(
        args.input_path, INPUT_FILE_FORMATS, check_readable=True
    )
    output_path:str = args.output_path
    common.validate_file_format(
        output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    common.print_start_message(args)

    genome, anno, *_ = common.load_references(args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    for record in parser.STARFusionParser.parse(fusion):
        if record.est_j < args.min_est_j:
            continue
        try:
            var_records = record.convert_to_variant_records(anno, genome)
        except err.GeneNotFoundError:
            continue
        variants.extend(var_records)

    logger.info('STAR-Fusion output %s loaded.', fusion)

    if not variants:
        logger.warning('No variant record is saved.')
        return

    genes_rank = anno.get_genes_rank()
    variants = sorted(variants, key=lambda x: genes_rank[x.location.seqname])

    logger.info('Variants sorted.')

    metadata = common.generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    logger.info('Variant info written to disk.')
