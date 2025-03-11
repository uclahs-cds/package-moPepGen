""" `parseRMATS` takes the alternative splicing event data called by
[rMATS](http://rnaseq-mats.sourceforge.net/) and converts them to a GVF file.
All five alternative splicing events are supported, including skipped exons,
alternative 5 splicing, alternative 3 splicing, mutually exclusive exons, and
retained introns. Both the tsv files with JC or JCEC suffix are supported.
The created GVF file can be then used to call for variant peptides using
[callVariant](call-variant.md)
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import argparse
from pathlib import Path
from moPepGen import get_logger, seqvar
from moPepGen.parser import RMATSParser
from moPepGen.cli import common


if TYPE_CHECKING:
    from typing import Dict, Set
    from logging import Logger

INPUT_FILE_FORMATS = ['.tsv', '.txt']
OUTPUT_FILE_FORMATS = ['.gvf']

def add_rmats_input_arg(parser:argparse.ArgumentParser, name:str, message:str,
        dest:str):
    """ add input arg for rMATS """
    message += f" Valid formats: {INPUT_FILE_FORMATS}"
    parser.add_argument(
        name,
        type=Path,
        help=message,
        metavar="<file>",
        default=None,
        dest=dest
    )


# pylint: disable=W0212
def add_subparser_parse_rmats(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseRMATs """
    ## parser_rmats
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='parseRMATS',
        help='Parse rMATS result for moPepGen to call variant peptides.',
        description='Parse the rMATS result to GVF format of variant'
        ' records for moPepGen to call variant peptides.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    add_rmats_input_arg(
        p, '--se',
        message="File path to the SE (skipped exons) junction count file output"
        " by rMATS. The file name should look like '*_SE.MATS.JC.txt' or "
        "'*_SE.MATS.JCEC.txt'.",
        dest='skipped_exon'
    )
    add_rmats_input_arg(
        p, '--a5ss',
        message="File path to the A5SS (alternative 5' splicint site) junction"
        " count file output by rMATS. The file name should look like"
        " '_S5SS.MATS.JC.txt' or '*_A5SS.MATS.JCEC.txt'.",
        dest='alternative_5_splicing'
    )
    add_rmats_input_arg(
        p, '--a3ss',
        message="File path to the A3SS (alternative 3' splicint site) junction"
        " count file output by rMATS. The file name should look like"
        " '_S3SS.MATS.JC.txt' or '*_A3SS.MATS.JCEC.txt'.",
        dest='alternative_3_splicing'
    )
    add_rmats_input_arg(
        p, '--mxe',
        message="File path to the MXE (mutually exclusive exons) junction"
        " count file output by rMATS. The file name should look like"
        " '_MXE.MATS.JC.txt' or '*_MXE.MATS.JCEC.txt'.",
        dest='mutually_exclusive_exons'
    )
    add_rmats_input_arg(
        p, '--ri',
        message="File path to the RI (retained intron) junction"
        " count file output by rMATS. The file name should look like"
        " '_RI.MATS.JC.txt' or '*_RI.MATS.JCEC.txt'.",
        dest='retained_intron'
    )
    p.add_argument(
        '--min-ijc',
        type=int,
        help='Minimal junction read count for the inclusion version to be'
        ' analyzed.',
        default=1
    )
    p.add_argument(
        '--min-sjc',
        type=int,
        help='Minimal junction read count for the skipped version to be'
        ' analyzed.',
        default=1
    )

    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_rmats)
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


def parse_rmats(args:argparse.Namespace) -> None:
    """ Parse rMATS results into TSV """
    logger = get_logger()
    tally = TallyTable(logger)

    skipped_exon = args.skipped_exon
    alternative_5 = args.alternative_5_splicing
    alternative_3 = args.alternative_3_splicing
    mutually_exclusive = args.mutually_exclusive_exons
    retained_intron = args.retained_intron
    output_path:Path = args.output_path

    for file in [skipped_exon, alternative_3, alternative_5, mutually_exclusive,
            retained_intron]:
        if file is not None:
            common.validate_file_format(
                file, INPUT_FILE_FORMATS, check_readable=True
            )
    common.validate_file_format(
        output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    common.print_start_message(args)

    genome, anno, *_ = common.load_references(args, load_canonical_peptides=False)

    variants:Dict[str,Set[seqvar.VariantRecord]] = {}
    rmats_outputs = [
        ('SE', skipped_exon), ('A5SS', alternative_5), ('A3SS', alternative_3),
        ('MXE', mutually_exclusive), ('RI', retained_intron)
    ]
    for event_type, path in rmats_outputs:
        if path:
            logger.info("Start parsing %s file %s", event_type, path)
            for record in RMATSParser.parse(path, event_type):
                tally.total += 1
                try:
                    var_records = record.convert_to_variant_records(
                        anno=anno, genome=genome,
                        min_ijc=args.min_ijc, min_sjc=args.min_sjc
                    )
                except:
                    logger.error(record.gene_id)
                    raise
                if var_records:
                    tally.succeed += 1
                else:
                    tally.skipped += 1
                for var_record in var_records:
                    tx_id = var_record.transcript_id
                    if tx_id not in variants:
                        variants[tx_id] = set()
                    variants[tx_id].add(var_record)

    if not variants:
        logger.warning('No variant record is saved.')
        return

    tx_rank = anno.get_transcript_rank()
    ordered_keys = sorted(variants.keys(), key=lambda x:tx_rank[x])
    variants_sorted = []
    for key in ordered_keys:
        val = list(variants[key])
        val.sort()
        variants_sorted.extend(val)

    logger.info('Variants sorted.')

    metadata = common.generate_metadata(args)
    seqvar.io.write(variants_sorted, output_path, metadata)

    logger.info('Variants written to disk.')

    tally.log()
