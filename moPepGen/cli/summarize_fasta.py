""" `summarizeFasta` takes a variant peptide FASTA file output by callVariant
and summarize the count of variant peptides of each source groups. This
summary can then guide the database splitting for tiered custom database
searching. """
from __future__ import annotations
import argparse
from contextlib import contextmanager
from pathlib import Path
import sys
from typing import IO
import matplotlib.pyplot as plt
from moPepGen.cli import common
from moPepGen.aa.PeptidePoolSummarizer import PeptidePoolSummarizer


GVF_FILE_FORMAT = ['.gvf']
FASTA_FILE_FORMAT = ['.fasta', '.fa']
OUTPUT_FILE_FORMATS = ['.txt', 'tsv']
OUTPUT_IMAGE_FORMATS = ['.pdf', '.jpg', '.jpeg', '.png']


# pylint: disable=W0212
def add_subparser_summarize_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen splitFasta """
    p:argparse.ArgumentParser = subparser.add_parser(
        name='summarizeFasta',
        help='Summarize the variant peptide calling results',
        description='Summarize the variant peptide calling results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '--gvf',
        type=Path,
        help='File path to GVF files. All GVF files must be generated by'
        f" moPepGen parsers. Valid formats: {GVF_FILE_FORMAT}",
        metavar='<files>',
        nargs='+'
    )
    p.add_argument(
        '--variant-peptides',
        type=Path,
        help='File path to the variant peptide FASTA database file. Must be'
        f" generated by moPepGen callVariant. Valid formats: {FASTA_FILE_FORMAT}",
        metavar='<file>'
    )
    p.add_argument(
        '--noncoding-peptides',
        type=Path,
        help='File path to the noncoding peptide FASTA database file. Must be'
        f" generated by moPepGen callNoncoding. Valid formats: {FASTA_FILE_FORMAT}",
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--alt-translation-peptides',
        type=Path,
        help='File path to the alt translation peptide FASTA file. Must be'
        f"generated by moPepGen callAltTranslation. Valid formats: {FASTA_FILE_FORMAT}"
    )
    p.add_argument(
        '--order-source',
        type=str,
        help='Order of sources, separate by comma. E.g., SNP,SNV,Fusion',
        metavar='<value>'
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help='File path to the output file. If not given, the summary table'
        f" is printed to stdout. Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--output-image',
        type=Path,
        help=f"File path to the output barplot. Valid formats: {OUTPUT_IMAGE_FORMATS}",
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--ignore-missing-source',
        action='store_true',
        help='Ignore the sources missing from input GVF.'
    )
    group_plot_scale = p.add_mutually_exclusive_group()
    group_plot_scale.add_argument(
        '--plot-normal-scale',
        action='store_true',
        help='Draw the summary bar plot in normal scale.'
    )
    group_plot_scale.add_argument(
        '--plot-log-scale',
        action='store_true',
        help='Draw the summary bar plot in log scale.'
    )

    common.add_args_cleavage(p, enzyme_only=True)
    common.add_args_reference(p, genome=False, proteome=True)
    common.add_args_quiet(p)
    common.print_help_if_missing_args(p)
    p.set_defaults(func=summarize_fasta)
    return p

@contextmanager
def output_context(file:Path) -> IO:
    """ Create context for summary output """
    if file is None:
        yield sys.stdout
    else:
        handle = open(file, 'wt')
        yield handle
        handle.close()

def summarize_fasta(args:argparse.Namespace) -> None:
    """ Summarize varaint peptide FASTA """
    for file in args.gvf:
        common.validate_file_format(
            file, GVF_FILE_FORMAT, check_readable=True
        )
    common.validate_file_format(
        args.variant_peptides, FASTA_FILE_FORMAT, check_readable=True
    )

    if args.output_path is not None:
        common.validate_file_format(
            args.output_path, OUTPUT_FILE_FORMATS, check_writable=True
        )
    if args.output_image is not None:
        common.validate_file_format(
            args.output_image, OUTPUT_IMAGE_FORMATS, check_writable=True
        )

    common.print_start_message(args)

    _, anno, *_ = common.load_references(
        args, load_genome=False, load_proteome=False,
        load_canonical_peptides=False, check_protein_coding=True
    )

    source_order = {val:i for i,val in  enumerate(args.order_source.split(','))}\
        if args.order_source else None

    summarizer = PeptidePoolSummarizer(
        order=source_order,
        ignore_missing_source=args.ignore_missing_source
    )

    for gvf in args.gvf:
        with open(gvf, 'rt') as handle:
            summarizer.update_label_map(handle)

    summarizer.append_order_internal_sources()

    with open(args.variant_peptides, 'rt') as handle:
        summarizer.load_database(handle)

    if args.noncoding_peptides:
        with open(args.noncoding_peptides, 'rt') as handle:
            summarizer.load_database(handle)

    if args.alt_translation_peptides:
        with open(args.alt_translation_peptides, 'rt') as handle:
            summarizer.load_database(handle)

    summarizer.count_peptide_source(anno, args.cleavage_rule)

    with output_context(args.output_path) as handle:
        summarizer.write_summary_table(handle)

    if args.output_image:
        if args.plot_log_scale:
            scale = 'log'
        elif args.plot_normal_scale:
            scale = 'normal'
        else:
            scale = None
        fig, ax = plt.subplots(figsize=(8, 8))
        summarizer.create_barplot(ax=ax, scale=scale)
        fig.savefig(args.output_image, bbox_inches="tight")
