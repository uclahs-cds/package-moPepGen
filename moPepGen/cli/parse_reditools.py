""" Module for REDItools parser """
from typing import Dict, List
import argparse
import pickle
import pathlib
from moPepGen import logger, gtf, seqvar, parser
from .common import add_args_reference, add_args_verbose, \
    print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_parse_reditools(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseREDItools """

    p = subparsers.add_parser(
        name='parseREDItools',
        help='Parse REDItools result for moPepGen to call variant peptides.',
        description='Parse the REDItools result to a TVF format of variant'
        'records for moPepGen to call variant peptides. The genome'
    )
    p.add_argument(
        '-t', '--reditools-table',
        type=str,
        help='Path to the REDItools output table.',
        metavar='',
        required=True
    )
    p.add_argument(
        '--transcript-id-column',
        type=int,
        help='The column index for transcript ID. If your REDItools table does'
        'not contains it, use the AnnotateTable.py from the REDItools'
        'package. Defaults to 16',
        default=16,
        metavar=''
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar='',
        required=True
    )
    add_args_reference(p)
    add_args_verbose(p)
    p.set_defaults(func=parse_reditools)
    print_help_if_missing_args(p)

def parse_reditools(args:argparse.Namespace) -> None:
    """ Parse REDItools output and save it in the TVF format. """
    # unpack args
    table_file = args.reditools_table
    transcript_id_column = args.transcript_id_column
    index_dir:str = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.tvf'
    verbose = args.verbose

    if verbose:
        logger('moPepGen parseREDItools started.')

    if index_dir:
        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            anno = pickle.load(handle)

        if verbose:
            logger('Index annotation loaded')

    else:
        annotation_gtf:str = args.annotation_gtf

        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

    variants:Dict[str, List[seqvar.VariantRecord]] = {}

    for record in parser.REDItoolsParser.parse(table_file, transcript_id_column):
        _vars = record.convert_to_variant_records(anno)
        for variant in _vars:
            transcript_id = variant.location.seqname
            if transcript_id not in variants:
                variants[transcript_id] = []
            variants[transcript_id].append(variant)

    if verbose:
        logger(f'REDItools table {table_file} loaded.')

    for records in variants.values():
        records.sort()

    if verbose:
        logger('Variants sorted.')

    if index_dir:
        reference_index = pathlib.Path(index_dir).absolute()
        annotation_gtf = None
    else:
        reference_index = None
        annotation_gtf = pathlib.Path(annotation_gtf).absolute()

    metadata = seqvar.TVFMetadata(
        parser='parseREDITools',
        reference_index=reference_index,
        genome_fasta=None,
        annotation_gtf=annotation_gtf
    )

    all_records = []
    for records in variants.values():
        all_records.extend(records)

    seqvar.io.write(all_records, output_path, metadata)

    if verbose:
        logger("Variants written to disk.")


if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.table_file = 'test/files/reditools/CPT0208690010_merged_chr22.txt'
    test_args.transcript_id_column = 16
    test_args.index_dir = 'test/files/downsampled_set/gencode_v36_index'
    test_args.output_prefix = 'test/files/reditools/CPT0208690010_merged_chr22.mop'
    test_args.verbose = True
    parse_reditools(test_args)
