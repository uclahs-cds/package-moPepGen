""" `parseREDItools` takes RNA editing results called by
[REDItools](https://github.com/BioinfoUNIBA/REDItools) and save as a GVF file.
The GVF file can then be used to call variant peptides using
[callVariant](call-variant.md)
"""
from __future__ import annotations
import argparse
from typing import Dict, List, TYPE_CHECKING
from moPepGen import logger, seqvar, parser
from .common import add_args_reference, add_args_verbose, add_args_source,\
    add_args_output_prefix, print_start_message,print_help_if_missing_args,\
    load_references, generate_metadata


# pylint: disable=W0212
def add_subparser_parse_reditools(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseREDItools """

    p = subparsers.add_parser(
        name='parseREDItools',
        help='Parse REDItools result for moPepGen to call variant peptides.',
        description='Parse the REDItools result to a GVF format of variant'
        'records for moPepGen to call variant peptides. The genome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-t', '--reditools-table',
        type=str,
        help='Path to the REDItools output table.',
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--transcript-id-column',
        type=int,
        help='The column index for transcript ID. If your REDItools table does'
        'not contains it, use the AnnotateTable.py from the REDItools'
        'package.',
        default=16,
        metavar='<number>'
    )
    add_args_output_prefix(p)
    add_args_source(p)
    add_args_reference(p, genome=False, proteome=False)
    add_args_verbose(p)
    p.set_defaults(func=parse_reditools)
    print_help_if_missing_args(p)
    return p

def parse_reditools(args:argparse.Namespace) -> None:
    """ Parse REDItools output and save it in the GVF format. """
    # unpack args
    table_file = args.reditools_table
    transcript_id_column = args.transcript_id_column
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    _, anno, *_ = load_references(args, load_genome=False, load_canonical_peptides=False)

    variants:Dict[str, List[seqvar.VariantRecord]] = {}

    for record in parser.REDItoolsParser.parse(table_file, transcript_id_column):
        _vars = record.convert_to_variant_records(anno)
        for variant in _vars:
            transcript_id = variant.location.seqname
            if transcript_id not in variants:
                variants[transcript_id] = []
            variants[transcript_id].append(variant)

    if args.verbose:
        logger(f'REDItools table {table_file} loaded.')

    for records in variants.values():
        records.sort()

    if args.verbose:
        logger('Variants sorted.')

    metadata = generate_metadata(args)

    all_records = []
    for records in variants.values():
        all_records.extend(records)

    seqvar.io.write(all_records, output_path, metadata)

    if args.verbose:
        logger("Variants written to disk.")


if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.table_file = 'test/files/reditools/CPT0208690010_merged_chr22.txt'
    test_args.transcript_id_column = 16
    test_args.index_dir = 'test/files/downsampled_set/gencode_v36_index'
    test_args.output_prefix = 'test/files/reditools/CPT0208690010_merged_chr22.mop'
    test_args.verbose = True
    parse_reditools(test_args)
