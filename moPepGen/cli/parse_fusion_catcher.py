""" `parseFusionCatcher` takes the identified fusion transcript results from
[FusionCatcher](https://github.com/ndaniel/fusioncatcher) and save as a
GVF file. The GVF file can be later used to call variant peptides using
[callVariant](call-variant.md)."""
from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path
import argparse
from moPepGen import get_logger, seqvar, parser, err
from moPepGen.cli import common


if TYPE_CHECKING:
    from typing import List
    from logging import Logger

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
    common.add_args_skip_failed(p)
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_fusion_catcher)
    common.print_help_if_missing_args(p)
    return p

class TallyTable():
    """ Tally table """
    def __init__(self, logger:Logger):
        """ Constructor """
        self.total:int = 0
        self.succeed:int = 0
        self.skipped:TallyTableSkipped = TallyTableSkipped()
        self.logger = logger

    def log(self):
        """ Show tally results """
        self.logger.info("Summary:")
        self.logger.info("Totally records read: %i", self.total)
        self.logger.info("Records successfully processed: %i", self.succeed)
        self.logger.info("Records skipped: %i", self.skipped.total)
        if self.skipped.total > 0:
            self.logger.info("Out of those skipped,")
            self.logger.info("    Invalid gene ID: %i", self.skipped.invalid_gene_id)
            self.logger.info("    Invalid position: %i", self.skipped.invalid_position)
            self.logger.info("    Insufficient evidence: %i", self.skipped.insufficient_evidence)

class TallyTableSkipped():
    """ Tally table for failed ones """
    def __init__(self):
        """ constructor """
        self.invalid_gene_id:int = 0
        self.invalid_position:int = 0
        self.insufficient_evidence:int = 0
        self.total:int = 0

def parse_fusion_catcher(args:argparse.Namespace) -> None:
    """ Parse FusionCatcher output and save it in GVF format. """
    logger = get_logger()
    tally = TallyTable(logger)
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
        tally.total += 1
        if record.counts_of_common_mapping_reads > args.max_common_mapping \
                or record.spanning_unique_reads < args.min_spanning_unique:
            tally.skipped.insufficient_evidence += 1
            tally.skipped.total += 1
            continue
        try:
            var_records = record.convert_to_variant_records(anno, genome)
            variants.extend(var_records)
            tally.succeed += 1
        except err.GeneNotFoundError:
            tally.skipped.invalid_gene_id += 1
            tally.skipped.total += 1
            continue
        except:
            if args.skip_failed:
                tally.skipped.total += 1
                tally.skipped.invalid_position += 1
                continue
            raise

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

    tally.log()
