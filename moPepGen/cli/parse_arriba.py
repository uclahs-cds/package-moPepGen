""" `parseArriba` takes the identified fusion transcript results from
[Arriba](https://github.com/suhrig/arriba) and saves as a GVF file. The GVF
file can be later used to call variant peptides using
[callVariant](call-variant.md)."""
from typing import List
from pathlib import Path
import argparse
from moPepGen import get_logger, seqvar, parser, err
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.tsv', '.txt']
OUTPUT_FILE_FORMATS = ['.gvf']

# pylint: disable=W0212
def add_subparser_parse_arriba(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseArriba """

    p:argparse.ArgumentParser = subparsers.add_parser(
        name='parseArriba',
        help='Parse Arriba result for moPepGen to call variant peptides.',
        description='Parse the Arriba result to GVF format of variant'
        ' records for moPepGen to call variant peptides.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=p, formats=INPUT_FILE_FORMATS,
        message="File path to Arriba's output TSV file."
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--min-split-read1',
        type=int,
        help='Minimal split_read1 value.',
        metavar='<value>',
        default=1
    )
    p.add_argument(
        '--min-split-read2',
        type=int,
        help='Minimal split_read2 value.',
        metavar='<value>',
        default=1
    )
    p.add_argument(
        '--min-confidence',
        type=str,
        choices=parser.ArribaParser.ArribaConfidence.levels.keys(),
        help='Minimal confidence value.',
        metavar='<choice>',
        default='medium'
    )
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_arriba)
    common.print_help_if_missing_args(p)
    return p

def parse_arriba(args:argparse.Namespace) -> None:
    """ Parse Arriba output and save it in GVF format. """
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

    min_split_read1:int = args.min_split_read1
    min_split_read2:int = args.min_split_read2
    min_confidence:str = args.min_confidence

    common.print_start_message(args)

    genome, anno, *_ = common.load_references(args=args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    with open(fusion, 'rt') as handle:
        for record in parser.ArribaParser.parse(handle):
            if not record.gene_id1 in anno.genes or not record.gene_id2 in anno.genes:
                continue
            if not record.is_valid(min_split_read1, min_split_read2, min_confidence):
                continue
            if record.transcript_on_antisense_strand(anno):
                continue
            try:
                var_records = record.convert_to_variant_records(anno, genome)
            except err.GeneNotFoundError:
                continue
            variants.extend(var_records)

    logger.info('Arriba output %s loaded.', fusion)

    if not variants:
        logger.warning('No variant record is saved.')
        return

    genes_rank = anno.get_genes_rank()
    variants = sorted(variants, key=lambda x: genes_rank[x.location.seqname])

    logger.info('Variants sorted.')

    metadata = common.generate_metadata(args)

    seqvar.io.write(variants, output_path, metadata)

    logger.info("Variants written to disk.")
