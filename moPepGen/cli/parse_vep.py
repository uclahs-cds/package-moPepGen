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
from moPepGen.err import TranscriptionStopSiteMutationError, \
    TranscriptionStartSiteMutationError, warning
from moPepGen import seqvar, logger
from moPepGen.cli.common import add_args_output_prefix, add_args_reference, \
    add_args_quiet, add_args_source, print_start_message, \
    print_help_if_missing_args, load_references, generate_metadata


# pylint: disable=W0212
def add_subparser_parse_vep(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseVEP """

    p:argparse.ArgumentParser = subparsers.add_parser(
        name='parseVEP',
        help='Parse VEP output for moPepGen to call variant peptides.',
        description="Parse VEP output tsv to the GVF format of variant records"
        "for moPepGen to call variant peptides. The genome assembly FASTA and"
        "annotation GTF must come from the same GENCODE/ENSEMBL version, and"
        "must the consistent with the VEP output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '-i', '--vep-txt',
        type=Path,
        nargs='+',
        help='Path to VEP result txt file.',
        metavar='<file>',
        required=True
    )
    add_args_output_prefix(p)
    add_args_source(p)
    add_args_reference(p, proteome=False)
    add_args_quiet(p)
    p.set_defaults(func=parse_vep)
    print_help_if_missing_args(p)
    return p

def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    # unpack args
    vep_files:List[Path] = args.vep_txt
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    genome, anno, *_ = load_references(args, load_canonical_peptides=False)

    vep_records:Dict[str, List[seqvar.VariantRecord]] = {}

    for vep_file in vep_files:
        opener = gzip.open if vep_file.suffix == '.gz' else open
        with opener(vep_file, 'rt') as handle:
            for record in VEPParser.parse(handle):
                transcript_id = record.feature

                if transcript_id not in vep_records.keys():
                    vep_records[transcript_id] = []

                try:
                    record = record.convert_to_variant_record(anno, genome)
                except TranscriptionStopSiteMutationError:
                    continue
                except TranscriptionStartSiteMutationError:
                    continue

                vep_records[transcript_id].append(record)

        if not args.quiet:
            logger(f'VEP file {vep_file} loaded.')

    if not vep_records:
        if not args.quiet:
            warning('No variant record is saved.')
        return

    for records in vep_records.values():
        records.sort()

    if not args.quiet:
        logger('VEP sorting done.')

    metadata = generate_metadata(args)

    tx_rank = anno.get_transcirpt_rank()
    ordered_keys = sorted(vep_records.keys(), key=lambda x:tx_rank[x])

    all_records = []
    for key in ordered_keys:
        all_records.extend(vep_records[key])

    seqvar.io.write(all_records, output_path, metadata)

    if not args.quiet:
        logger('Variant info written to disk.')
