""" Module for CIRCexplorer parser """
from __future__ import annotations
from typing import List, Dict, TYPE_CHECKING
from pathlib import Path
from moPepGen import logger, circ
from moPepGen.parser import CIRCexplorerParser
from .common import add_args_reference, add_args_verbose, add_args_source,\
    print_start_message,print_help_if_missing_args, load_references, \
    generate_metadata


if TYPE_CHECKING:
    import argparse

# pylint: disable=W0212
def add_subparser_parse_circexplorer(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseCIRCexplorer """
    p = subparsers.add_parser(
        name='parseCIRCexplorer',
        help='Parse CIRCexplorer result',
        description='Parse CIRCexplorer result to a TSV format for moPepGen to'
        ' call variant peptides'
    )
    p.add_argument(
        '-i', '--input-path',
        type=Path,
        help='The input file path for CIRCexplorer result. Only the known'
        'circRNA result is supported.',
        required=True,
        metavar=''
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Output prefix',
        required=True,
        metavar=''
    )
    p.add_argument(
        '--min-read-number',
        type=int,
        help='Minimal number of junction read counts. Defaults to 2',
        default=2,
        metavar=''
    )
    add_args_source(p)
    add_args_reference(p, genome=False, proteome=False)
    add_args_verbose(p)
    p.set_defaults(func=parse_circexplorer)
    print_help_if_missing_args(p)

def parse_circexplorer(args:argparse.Namespace):
    """ Parse circexplorer known circRNA results. """
    input_path = args.input_path
    output_prefix = args.output_prefix
    output_path = output_prefix + '.tsv'

    print_start_message(args)

    _, anno, *_ = load_references(args, False, False)

    circ_records:Dict[str, List[circ.CircRNAModel]] = {}

    for record in CIRCexplorerParser.parse(input_path):
        if record.read_number < args.min_read_number:
            continue
        circ_record = record.convert_to_circ_rna(anno)
        gene_id = circ_record.gene_id
        if gene_id not in circ_records:
            circ_records[gene_id] = []
        circ_records[gene_id].append(circ_record)

    records = []
    for val in circ_records.values():
        records.extend(val)

    metadata = generate_metadata(args)

    circ.io.write(records, metadata, output_path)

    if args.verbose:
        logger("CircRNA records written to disk.")
