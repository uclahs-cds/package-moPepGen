""" VEP2VariantPeptides module """
from __future__ import annotations
from typing import Dict, List, TYPE_CHECKING
from pathlib import Path
from moPepGen.parser import VEPParser
from moPepGen.err import TranscriptionStopSiteMutationError
from moPepGen import seqvar, logger
from moPepGen.cli.common import add_args_reference, add_args_verbose, add_args_source,\
    print_start_message, print_help_if_missing_args, load_references, \
    generate_metadata


if TYPE_CHECKING:
    import argparse

# pylint: disable=W0212
def add_subparser_parse_vep(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseVEP """

    p = subparsers.add_parser(
        name='parseVEP',
        help='Parse VEP output for moPepGen to call variant peptides.',
        description="Parse VEP output tsv to the GVF format of variant records"
        "for moPepGen to call variant peptides. The genome assembly FASTA and"
        "annotation GTF must come from the same GENCODE/ENSEMBL version, and"
        "must the consistent with the VEP output."
    )

    p.add_argument(
        '-i', '--vep-txt',
        type=Path,
        nargs='+',
        help='Path to VEP result txt file.',
        metavar='',
        required=True
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar='',
        required=True
    )
    add_args_source(p)
    add_args_reference(p, proteome=False)
    add_args_verbose(p)
    p.set_defaults(func=parse_vep)
    print_help_if_missing_args(p)

def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    # unpack args
    vep_files:List[str] = args.vep_txt
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    genome, anno, *_ = load_references(args, load_canonical_peptides=False)

    vep_records:Dict[str, List[seqvar.VariantRecord]] = {}

    for vep_file in vep_files:
        for record in VEPParser.parse(vep_file):
            transcript_id = record.feature

            if transcript_id not in vep_records.keys():
                vep_records[transcript_id] = []

            try:
                record = record.convert_to_variant_record(anno, genome)
            except TranscriptionStopSiteMutationError:
                continue

            vep_records[transcript_id].append(record)

        if args.verbose:
            logger(f'VEP file {vep_file} loaded.')

    for records in vep_records.values():
        records.sort()

    if args.verbose:
        logger('VEP sorting done.')

    metadata = generate_metadata(args)

    all_records = []
    for records in vep_records.values():
        all_records.extend(records)

    seqvar.io.write(all_records, output_path, metadata)

    if args.verbose:
        logger('Variant info written to disk.')


if __name__ == '__main__':
    import argparse
    test_args = argparse.Namespace()
    test_args.command = 'parseVEP'
    test_args.vep_txt = [
        'test/files/CPCG0100_gencode_indel_ENST00000516353.1.txt'
    ]
    test_args.index_dir = 'test/files/gencode_34_index'
    test_args.genome_fasta = None
    test_args.annotation_gtf = None
    test_args.output_prefix = 'test/files/vep/vep'
    test_args.verbose = True
    parse_vep(args=test_args)
