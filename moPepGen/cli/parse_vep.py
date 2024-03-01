""" `parseVEP` takes the output of Ensembl's [Variant Effector Predictor
](https://uswest.ensembl.org/info/docs/tools/vep/index.html) (VEP) and
convert it into the GVF file format that moPepGen internally uses. The result
VEP file can then be parsed to moPepGen's `callVariant` subcommand to call for
variant peptide sequences.
"""
from __future__ import annotations
import argparse
import gzip
from typing import Dict, List
from pathlib import Path
from moPepGen.parser import VEPParser
from moPepGen.err import MNVParsingError, TranscriptionStopSiteMutationError, \
    TranscriptionStartSiteMutationError, warning
from moPepGen import seqvar, get_logger
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.tsv', '.txt', '.tsv.gz', '.txt.gz']
OUTPUT_FILE_FORMATS = ['.gvf']

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
        message='File path to the VEP output TXT file.'
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    common.add_args_source(p)
    common.add_args_reference(p, proteome=False)
    common.add_args_debug_level(p)
    p.set_defaults(func=parse_vep)
    common.print_help_if_missing_args(p)
    return p

def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    logger = get_logger()
    # unpack args
    vep_files:List[Path] = args.input_path
    for file in vep_files:
        common.validate_file_format(
            file, INPUT_FILE_FORMATS, check_readable=True
        )
    output_path:Path = args.output_path
    common.validate_file_format(
        output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    common.print_start_message(args)

    genome, anno, *_ = common.load_references(args, load_canonical_peptides=False)

    vep_records:Dict[str, List[seqvar.VariantRecord]] = {}

    for vep_file in vep_files:
        opener = gzip.open if vep_file.suffix == '.gz' else open
        with opener(vep_file, 'rt') as handle:
            for record in VEPParser.parse(handle):
                transcript_id = record.feature

                if transcript_id not in vep_records:
                    vep_records[transcript_id] = []

                try:
                    record = record.convert_to_variant_record(anno, genome)
                except TranscriptionStopSiteMutationError:
                    continue
                except TranscriptionStartSiteMutationError:
                    continue
                except MNVParsingError:
                    warning(
                        f"MNVs are not currently supported. Skipping record: {record}"
                    )
                    continue

                vep_records[transcript_id].append(record)

        logger.info('VEP file %s loaded.', vep_file)

    if not vep_records:
        logger.warning('No variant record is saved.')
        return

    for records in vep_records.values():
        records.sort()

    logger.info('VEP sorting done.')

    metadata = common.generate_metadata(args)

    tx_rank = anno.get_transcript_rank()
    ordered_keys = sorted(vep_records.keys(), key=lambda x:tx_rank[x])

    all_records = []
    for key in ordered_keys:
        all_records.extend(vep_records[key])

    seqvar.io.write(all_records, output_path, metadata)

    logger.info('Variant info written to disk.')
