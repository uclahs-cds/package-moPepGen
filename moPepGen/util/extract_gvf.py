""" Extract records from GVF files """
import argparse
from pathlib import Path
from typing import List
from moPepGen import seqvar
from moPepGen import circ
from moPepGen.cli.common import print_help_if_missing_args, add_args_debug_level
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen.circ import CircRNAModel


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ add args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='extractGVF',
        help='Extract records from GVF files'
    )
    parser.add_argument(
        '-i', '--input-gvf',
        type=Path,
        help='Input GVF files',
        nargs='*',
        metavar='<files>',
        required=True
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Output directory',
        metavar='<file>',
        required=True
    )
    parser.add_argument(
        '--gene-list',
        type=str,
        help='Gene list to filter.',
        nargs='*',
        metavar='<values>',
        default=[]
    )
    parser.add_argument(
        '--tx-list',
        type=str,
        help='Transcript list to filter.',
        nargs='*',
        metavar='<values>',
        default=[]
    )
    parser.set_defaults(func=main)
    add_args_debug_level(parser)
    print_help_if_missing_args(parser)
    return parser

def is_wanted_variant(variant:seqvar.VariantRecord, tx_list:List[str],
        gene_list:List[str]) -> bool:
    """ check if the variant record is wanted """
    return variant.location.seqname in gene_list or \
        variant.attrs['TRANSCRIPT_ID'] in tx_list

def is_wanted_circ_rna(variant:CircRNAModel, tx_list:List[str],
        gene_list:List[str]) -> bool:
    """ check if the circRNA record is wanted """
    return variant.gene_id in gene_list or variant.transcript_id in tx_list

def main(args:argparse.Namespace):
    """ Extract records from GVF files """
    gvf_files:List[Path] = args.input_gvf
    gene_list:List[str] = args.gene_list
    tx_list:List[str] = args.tx_list
    output_dir:Path = args.output_dir

    for gvf_file in gvf_files:
        with open(gvf_file, 'rt') as handle:
            metadata = GVFMetadata.parse(handle)
            output_file = output_dir/f"{metadata.source}.gvf"
            if metadata.is_circ_rna():
                records:List[CircRNAModel] = []
                for record in circ.io.parse(handle):
                    if is_wanted_circ_rna(record, tx_list, gene_list):
                        records.append(record)
                with open(output_file, 'wt') as out_handle:
                    circ.io.write(records, metadata, out_handle)
            else:
                records:List[seqvar.VariantRecord] = []
                for record in seqvar.io.parse(handle):
                    if is_wanted_variant(record, tx_list, gene_list):
                        records.append(record)
                seqvar.io.write(records, output_file, metadata)
