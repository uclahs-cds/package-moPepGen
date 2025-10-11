""" `parseVEP` takes the output of Ensembl's [Variant Effector Predictor
](https://uswest.ensembl.org/info/docs/tools/vep/index.html) (VEP) and
convert it into the GVF file format that moPepGen internally uses. The result
VEP file can then be parsed to moPepGen's `callVariant` subcommand to call for
variant peptide sequences.
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import argparse
import gzip
from pathlib import Path
from moPepGen.parser import VEPParser
from moPepGen.err import TranscriptionStopSiteMutationError, TranscriptionStartSiteMutationError
from moPepGen import seqvar, get_logger
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.tsv', '.txt', '.tsv.gz', '.txt.gz', '.vcf', '.vcf.gz']
OUTPUT_FILE_FORMATS = ['.gvf']

if TYPE_CHECKING:
    from typing import Dict, List, Tuple, Set
    from logging import Logger

# pylint: disable=W0212
def add_subparser_parse_vep(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseVEP """

    p:argparse.ArgumentParser = subparsers.add_parser(
        name='parseVEP',
        help='Parse VEP output for moPepGen to call variant peptides.',
        description='Parse VEP output tsv to the GVF format of variant records'
        ' for moPepGen to call variant peptides. The genome assembly FASTA and'
        ' annotation GTF must come from the same GENCODE/ENSEMBL version, and'
        ' must the consistent with the VEP output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=p, formats=INPUT_FILE_FORMATS, plural=True,
        message='File path to the VEP output TXT or VCF file.'
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help=f'File path to the output file. Valid formats: {OUTPUT_FILE_FORMATS}',
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--output-prefix',
        type=Path,
        help='Output prefix. Only used when inputs are VCF files.',
        default=None,
        metavar = '<value>'
    )
    p.add_argument(
        '--samples',
        type=str,
        nargs='+',
        default=[],
        help='Samples to be parsed from the VCF file. If not provided, all samples'
        ' from the VCF file will be parsed. This option is only used when the inputs'
        ' are VCF files.'
    )
    common.add_args_source(p)
    common.add_args_skip_failed(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_vep)
    common.print_help_if_missing_args(p)
    return p

class TallyTable():
    """ Tally table """
    def __init__(self, logger:Logger):
        """ Constructor """
        self.total:int = 0
        self.succeed:int = 0
        self.failed:TallyTableFailed = TallyTableFailed()
        self.logger = logger

    def log(self):
        """ Show tally results """
        self.logger.info("Totally records read: %i", self.total)
        self.logger.info("Records successfully processed: %i", self.succeed)
        self.logger.info("Records failed: %i", self.failed.total)
        if self.failed.total > 0:
            self.logger.info("Out of those failed,")
            self.logger.info("Start codon mutation: %i", self.failed.start_site_mutation)
            self.logger.info("Stop codon mutation: %i", self.failed.stop_site_mutation)

class TallyTableFailed():
    """ Tally table for failed ones """
    def __init__(self):
        """ constructor """
        self.start_site_mutation:int = 0
        self.stop_site_mutation:int = 0
        self.total:int = 0

def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    logger = get_logger()
    # unpack args
    vep_files:List[Path] = args.input_path
    format = None
    for file in vep_files:
        common.validate_file_format(
            file, INPUT_FILE_FORMATS, check_readable=True
        )
        if '.vcf' in file.suffixes:
            format_i = 'vcf'
        else:
            format_i = 'tsv'
        if format is None:
            format = format_i
        elif format != format_i:
            raise ValueError(
                "All input files must be in the same format. "
                f"Got {format} and {format_i}."
            )
    if format == 'tsv':
        output_path:Path = args.output_path
        if output_path is None:
            raise ValueError(
                "--output-path must be provided when the input files are TSV or TXT."
            )
        common.validate_file_format(
            output_path, OUTPUT_FILE_FORMATS, check_writable=True
        )
    else:
        output_prefix:Path = args.output_prefix
        if output_prefix is None:
            raise ValueError(
                "--output-prefix must be provided when the input files are VCF."
            )
        output_dir = output_prefix.parent
        output_prefix = output_prefix.name

    common.print_start_message(args)

    ref_data = common.load_references(args, load_canonical_peptides=False)
    genome = ref_data.genome
    anno = ref_data.anno

    # sample -> transcript_id -> list of variant records
    vep_records:Dict[str, Dict[str, List[seqvar.VariantRecord]]] = {}

    if format == 'vcf':
        samples = args.samples
    else:
        samples = ['default']

    tally = TallyTable(logger)

    phase_sets:Set[Tuple[str, str]] = set()

    for vep_file in vep_files:
        opener = gzip.open if vep_file.suffix == '.gz' else open
        with opener(vep_file, 'rt') as handle:
            it = VEPParser.parse(
                handle=handle,
                format=format,
                samples=samples,
                current_phase_sets=phase_sets
            )
            for vep_record in it:
                tally.total += 1
                sample = vep_record.extra.get('SAMPLE', 'default')
                transcript_id = vep_record.feature

                try:
                    record = vep_record.convert_to_variant_record(anno, genome)
                    tally.succeed += 1
                except TranscriptionStopSiteMutationError:
                    tally.failed.total += 1
                    tally.failed.stop_site_mutation += 1
                    continue
                except TranscriptionStartSiteMutationError:
                    tally.failed.total += 1
                    tally.failed.start_site_mutation += 1
                    continue
                except:
                    if args.skip_failed:
                        logger.warning(
                            "VEP record failed to convert: %s", str(vep_record)
                        )
                        tally.failed.total += 1
                        continue
                    raise

                if sample not in vep_records:
                    vep_records[sample] = {}
                if transcript_id not in vep_records[sample]:
                    vep_records[sample][transcript_id] = []
                vep_records[sample][transcript_id].append(record)

        logger.info('VEP file %s loaded.', vep_file)

    tally.log()

    if not vep_records:
        logger.warning('No variant record is saved.')
        return

    for sample_records in vep_records.values():
        for records in sample_records.values():
            records.sort()

    logger.info('GVF records sorted.')

    tx_rank = anno.get_transcript_rank()

    for sample, sample_records in vep_records.items():
        ordered_keys = sorted(sample_records.keys(), key=lambda x: tx_rank[x])

        sorted_records = []
        for key in ordered_keys:
            sorted_records.extend(sample_records[key])

        if format == 'vcf':
            output_path = output_dir/f"{output_prefix}_{sample}.gvf"

        metadata = common.generate_metadata(args, phase_pairs=phase_sets)
        seqvar.io.write(sorted_records, output_path, metadata)

        logger.info('Variant GVF written: %s', output_path)
