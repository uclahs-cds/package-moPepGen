""" `parseFusionCatcher` takes the identified fusion transcript results from
[FusionCatcher](https://github.com/ndaniel/fusioncatcher) and save as a
GVF file. The GVF file can be later used to call variant peptides using
[callVariant](call-variant.md)."""
from typing import List
from pathlib import Path
import argparse
from moPepGen import get_logger, seqvar, parser, err
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.tsv', '.txt']
OUTPUT_FILE_FORMATS = ['.gvf']

# pylint: disable=W0212
def add_subparser_parse_fusion_catcher(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseFusionCatcher """

    p = subparsers.add_parser(
        name='parseFusionCatcher',
        help='Parse FusionCatcher result for moPepGen to call variant peptides.',
        description='Parse the FusionCatcher result to GVF format of variant'
        ' records for moPepGen to call variant peptides. The genome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=p, formats=INPUT_FILE_FORMATS,
        message="File path to FusionCatcher's output TSV file."
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--max-common-mapping',
        type=int,
        help='Maximal number of common mapping reads.',
        metavar='<number>',
        default=0
    )
    p.add_argument(
        '--min-spanning-unique',
        help='Minimal spanning unique reads.',
        type=int,
        default=5,
        metavar='<number>'
    )
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_fusion_catcher)
    common.print_help_if_missing_args(p)
    return p

def parse_fusion_catcher(args:argparse.Namespace) -> None:
    """ Parse FusionCatcher output and save it in GVF format. """
    logger = get_logger()
    # unpack args
    fusion = args.input_path
    output_path:Path = args.output_path
    common.validate_file_format(
        fusion, INPUT_FILE_FORMATS, check_readable=True
    )
    common.validate_file_format(
        output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    common.print_start_message(args)

    genome, anno, *_ = common.load_references(args=args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    for record in parser.FusionCatcherParser.parse(fusion):
        if record.counts_of_common_mapping_reads > args.max_common_mapping:
            continue
        if record.spanning_unique_reads < args.min_spanning_unique:
            continue
        try:
            var_records = record.convert_to_variant_records(anno, genome)
        except err.GeneNotFoundError:
            continue
        variants.extend(var_records)

    logger.info('FusionCatcher output %s loaded.', fusion)

    if not variants:
        logger.warning('No variant record is saved.')
        return

    genes_rank = anno.get_genes_rank()
    variants = sorted(variants, key=lambda x: genes_rank[x.location.seqname])

    logger.info('Variants sorted.')

    metadata = common.generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    logger.info("Variants written to disk.")
