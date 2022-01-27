""" `parseRMATS` takes the alternative splicing event data called by
[rMATS](http://rnaseq-mats.sourceforge.net/) and converts them to a GVF file.
All five alternative splicing events are supported, including skipped exons,
alternative 5 splicing, alternative 3 splicing, mutually exclusive exons, and
retained introns. Both the tsv files with JC or JCEC suffix are supported.
The created GVF file can be then used to call for variant peptides using
[callVariant](call-variant.md)
"""
from __future__ import annotations
import argparse
from logging import warning
from typing import Dict, Set
from pathlib import Path
from moPepGen import logger, seqvar
from moPepGen.parser import RMATSParser
from .common import add_args_output_prefix, add_args_reference, \
    add_args_quiet, add_args_source, print_start_message, \
    print_help_if_missing_args, load_references, generate_metadata


# pylint: disable=W0212
def add_subparser_parse_rmats(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseRMATs """
    ## parser_rmats
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='parseRMATS',
        help='Parse rMATS result for moPepGen to call variant peptides.',
        description='Parse the rMATS result to GVF format of variant'
        'records for moPepGen to call variant peptides.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '--se',
        type=Path,
        help="Skipped exon junction count txt file. The file name should have"
        " the pattern of *_SE.MATS.JC.txt or *_SE.MATS.JCEC.txt",
        metavar='<file>',
        default=None,
        dest='skipped_exon'
    )
    p.add_argument(
        '--a5ss',
        type=Path,
        help="Alternative 5' splicing junction count txt file. The file name"
        ' should have the pattern of *_S5SS.MATS.JC.txt or *_A5SS.MATS.JCEC.txt',
        metavar='<file>',
        default=None,
        dest='alternative_5_splicing'
    )
    p.add_argument(
        '--a3ss',
        type=Path,
        help="Alternative 3' splicing junction count txt file. The file name"
        ' should have the pattern of *_S3SS.MATS.JC.txt or *_A3SS.MATS.JCEC.txt',
        metavar='<file>',
        default=None,
        dest='alternative_3_splicing'
    )
    p.add_argument(
        '--mxe',
        type=Path,
        help="Mutually exclusive junction count txt file. The file name should"
        " have the pattern of *_MXE.MATS.JC.txt or *_MXE.MATS.JCEC.txt',",
        metavar='<file>',
        default=None,
        dest='mutually_exclusive_exons'
    )
    p.add_argument(
        '--ri',
        type=Path,
        help="Retained intron junction count txt file. The file name should"
        " have the pattern of *_RI.MATS.JC.txt or *_RI.MATS.JCEC.txt",
        metavar='<file>',
        default=None,
        dest='retained_intron'
    )
    p.add_argument(
        '--min-ijc',
        type=int,
        help="Minimal junction read count for the inclusion version to be"
        " analyzed.",
        default=1
    )
    p.add_argument(
        '--min-sjc',
        type=int,
        help="Minimal junction read count for the skipped version to be"
        " analyzed.",
        default=1
    )

    add_args_output_prefix(p)
    add_args_source(p)
    add_args_reference(p, proteome=False)
    add_args_quiet(p)
    p.set_defaults(func=parse_rmats)
    print_help_if_missing_args(p)
    return p

def parse_rmats(args:argparse.Namespace) -> None:
    """ Parse rMATS results into TSV """
    skipped_exon = args.skipped_exon
    alternative_5 = args.alternative_5_splicing
    alternative_3 = args.alternative_3_splicing
    mutually_exclusive = args.mutually_exclusive_exons
    retained_intron = args.retained_intron
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.gvf'

    print_start_message(args)

    genome, anno, *_ = load_references(args, load_canonical_peptides=False)

    variants:Dict[str,Set[seqvar.VariantRecord]] = {}
    rmats_outputs = [
        ('SE', skipped_exon), ('A5SS', alternative_5), ('A3SS', alternative_3),
        ('MXE', mutually_exclusive), ('RI', retained_intron)
    ]
    for event_type, path in rmats_outputs:
        if path:
            logger(f"Start parsing {event_type} file {path}")
            for record in RMATSParser.parse(path, event_type):
                try:
                    var_records = record.convert_to_variant_records(
                        anno=anno, genome=genome,
                        min_ijc=args.min_ijc, min_sjc=args.min_sjc
                    )
                except:
                    logger(record.gene_id)
                    raise
                for var_record in var_records:
                    gene_id = var_record.location.seqname
                    if gene_id not in variants:
                        variants[gene_id] = set()
                    variants[gene_id].add(var_record)

    if not variants:
        if not args.quiet:
            warning('No variant record is saved.')
        return

    genes_rank = anno.get_genes_rank()
    ordered_keys = sorted(variants.keys(), key=lambda x:genes_rank[x])
    variants_sorted = []
    for key in ordered_keys:
        val = list(variants[key])
        val.sort()
        variants_sorted.extend(val)

    if not args.quiet:
        logger('Variants sorted.')

    metadata = generate_metadata(args)
    seqvar.io.write(variants_sorted, output_path, metadata)

    if not args.quiet:
        logger("Variants written to disk.")
