""" `parseREDItools` takes RNA editing results called by
[REDItools](https://github.com/BioinfoUNIBA/REDItools) and saves them as a GVF
file. The GVF file can then be used to call variant peptides using
[callVariant](call-variant.md)
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import argparse
from pathlib import Path
from moPepGen import get_logger, seqvar, parser
from moPepGen.cli import common


if TYPE_CHECKING:
    from typing import Dict, List
    from logging import Logger

INPUT_FILE_FORMATS = ['.tsv', '.txt']
OUTPUT_FILE_FORMATS = ['.gvf']

# pylint: disable=W0212
def add_subparser_parse_reditools(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseREDItools """

    p:argparse.ArgumentParser = subparsers.add_parser(
        name='parseREDItools',
        help='Parse REDItools result for moPepGen to call variant peptides.',
        description='Parse the REDItools result to a GVF format of variant'
        ' records for moPepGen to call variant peptides. The genome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=p, formats=INPUT_FILE_FORMATS,
        message="File path to REDItools' TSV output."
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--transcript-id-column',
        type=int,
        help='The column index for transcript ID. If your REDItools table does'
        ' not contains it, use the AnnotateTable.py from the REDItools'
        ' package.',
        default=17,
        metavar='<number>'
    )
    p.add_argument(
        '--min-coverage-alt',
        type=int,
        help='Minimal read coverage of alterations to be parsed.',
        default=3,
        metavar='<number>'
    )
    p.add_argument(
        '--min-frequency-alt',
        type=float,
        help='Minimal frequency of alteration to be parsed.',
        default=0.1,
        metavar='<value>'
    )
    p.add_argument(
        '--min-coverage-rna',
        type=int,
        help='Minimal read coverage at the alteration site of RNAseq data of'
        ' reference and all alterations.',
        default=10,
        metavar='<value>'
    )
    p.add_argument(
        '--min-coverage-dna',
        type=int,
        help='Minimal read coverage at the alteration site of WGS. Set it to'
        ' -1 to skip checking this.',
        default=10,
        metavar='<number>'
    )
    common.add_args_source(p)
    common.add_args_reference(p, genome=False, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_reditools)
    common.print_help_if_missing_args(p)
    return p

class TallyTable():
    """ Tally table """
    def __init__(self, logger:Logger):
        """ Constructor """
        self.total:int = 0
        self.succeed:int = 0
        self.skipped:int = 0
        self.logger = logger

    def log(self):
        """ Show tally results """
        self.logger.info("Summary:")
        self.logger.info("Totally records read: %i", self.total)
        self.logger.info("Records successfully processed: %i", self.succeed)
        self.logger.info("Records skipped: %i", self.skipped)

def parse_reditools(args:argparse.Namespace) -> None:
    """ Parse REDItools output and save it in the GVF format. """
    logger = get_logger()
    tally = TallyTable(logger)
    # unpack args
    table_file:Path = args.input_path
    output_path:Path = args.output_path
    common.validate_file_format(
        table_file, INPUT_FILE_FORMATS, check_readable=True
    )
    common.validate_file_format(
        output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    transcript_id_column = args.transcript_id_column - 1
    min_coverage_alt:int = args.min_coverage_alt
    min_frequency_alt:int = args.min_frequency_alt
    min_coverage_rna:int = args.min_coverage_rna
    min_coverage_dna:int = args.min_coverage_dna

    common.print_start_message(args)

    _, anno, *_ = common.load_references(args, load_genome=False, load_canonical_peptides=False)

    variants:Dict[str, List[seqvar.VariantRecord]] = {}

    for record in parser.REDItoolsParser.parse(table_file, transcript_id_column):
        tally.total += 1
        _vars = record.convert_to_variant_records(
            anno=anno,
            min_coverage_alt=min_coverage_alt,
            min_frequency_alt=min_frequency_alt,
            min_coverage_rna=min_coverage_rna,
            min_coverage_dna=min_coverage_dna
        )
        if not _vars:
            tally.skipped += 1
        else:
            tally.succeed += 1
        for variant in _vars:
            gene_id = variant.location.seqname
            if gene_id not in variants:
                variants[gene_id] = []
            variants[gene_id].append(variant)

    logger.info('REDItools table %s loaded.', table_file)

    if not variants:
        logger.warning('No variant record is saved.')
        return

    for records in variants.values():
        records.sort()

    logger.info('Variants sorted.')

    metadata = common.generate_metadata(args)

    genes_rank = anno.get_genes_rank()
    ordered_keys = sorted(variants.keys(), key=lambda x:genes_rank[x])

    all_records = []
    for key in ordered_keys:
        records = variants[key]
        all_records.extend(records)

    seqvar.io.write(all_records, output_path, metadata)

    logger.info('Variants written to disk.')

    tally.log()
