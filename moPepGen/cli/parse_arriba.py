""" `parseArriba` takes the identified fusion transcript results from
[Arriba](https://github.com/suhrig/arriba) and saves as a GVF file. The GVF
file can be later used to call variant peptides using
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
    common.add_args_skip_failed(p)
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_arriba)
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
            self.logger.info("    Antisense strand: %i", self.skipped.antisense_strand)

class TallyTableSkipped():
    """ Tally table for failed ones """
    def __init__(self):
        """ constructor """
        self.invalid_gene_id:int = 0
        self.invalid_position:int = 0
        self.insufficient_evidence:int = 0
        self.antisense_strand:int = 0
        self.total:int = 0

def parse_arriba(args:argparse.Namespace) -> None:
    """ Parse Arriba output and save it in GVF format. """
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

    min_split_read1:int = args.min_split_read1
    min_split_read2:int = args.min_split_read2
    min_confidence:str = args.min_confidence

    common.print_start_message(args)

    genome, anno, *_ = common.load_references(args=args, load_canonical_peptides=False)

    variants:List[seqvar.VariantRecord] = []

    with open(fusion, 'rt') as handle:
        for record in parser.ArribaParser.parse(handle):
            tally.total += 1
            if not record.gene_id1 in anno.genes or not record.gene_id2 in anno.genes:
                tally.skipped.invalid_gene_id += 1
                tally.skipped.total += 1
                continue
            if not record.is_valid(min_split_read1, min_split_read2, min_confidence):
                tally.skipped.insufficient_evidence += 1
                tally.skipped.total += 1
                continue
            if record.transcript_on_antisense_strand(anno):
                tally.skipped.antisense_strand += 1
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
                    tally.skipped.invalid_position += 1
                    tally.skipped.total += 1
                    continue
                raise

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

    tally.log()
