""" `splitFasta` takes the FASTA file with variant peptide sequences called
by [`callVariant`](./call-variant.md) with or without noncoding novel peptides
called by [`callNoncoding`](./call-noncoding.md), and splits peptide sequences
into databases. The split database FASTA files can be used for sequential
library searching. """
from __future__ import annotations
import argparse
from pathlib import Path
from moPepGen.aa import PeptidePoolSplitter
from moPepGen import SPLIT_DATABASE_KEY_SEPARATER, logger
from moPepGen.cli import common


GVF_FILE_FORMAT = ['.gvf']
FASTA_FILE_FORMAT = ['.fasta', '.fa']


# pylint: disable=W0212
def add_subparser_split_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen splitFasta """
    p:argparse.ArgumentParser = subparser.add_parser(
        name='splitFasta',
        help='Split variant peptide FASTA database generated by moPepGen.',
        description='Split variant peptide FASTA database generated by'
        ' moPepGen into separate files.',
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
        f"generated by moPepGen callAltTranslation. Valid formats: {FASTA_FILE_FORMAT}",
        metavar='<file>'
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=Path,
        help='Output prefix',
        metavar = '<value>'
    )
    p.add_argument(
        '--order-source',
        type=str,
        help='Order of sources, separate by comma. E.g., SNP,SNV,Fusion',
        metavar='<value>'
    )
    p.add_argument(
        '--group-source',
        type=str,
        help='Group sources. The peptides with sources grouped will be written'
        ' to the same FASTA file. E.g., "PointMutation:gSNP,sSNV'
        ' INDEL:gINDEL,sINDEL".',
        metavar='<value>',
        nargs='*'
    )
    p.add_argument(
        '--max-source-groups',
        type=int,
        help='Maximal number of different source groups to be separate into'
        ' individual database FASTA files. Defaults to 1',
        default=1,
        metavar='<number>'
    )
    p.add_argument(
        '--additional-split',
        type=str,
        help='For peptides that were not already split into FASTAs up to'
        ' max_source_groups, those involving the following source will be split'
        ' into additional FASTAs with decreasing priority. E.g., '
        " 'gSNP-Noncoding', 'gSNP-Noncoding gSNP-gINDEL'",
        metavar='<value>',
        default=None,
        nargs="*"
    )

    common.add_args_reference(p, genome=False, proteome=True)
    common.add_args_quiet(p)
    common.print_help_if_missing_args(p)
    p.set_defaults(func=split_fasta)
    return p

def split_fasta(args:argparse.Namespace) -> None:
    """ Split peptide database """
    for file in args.gvf:
        common.validate_file_format(
            file, GVF_FILE_FORMAT, check_readable=True
        )
    for file in [args.variant_peptides, args.noncoding_peptides]:
        if file is not None:
            common.validate_file_format(
                file, FASTA_FILE_FORMAT, check_readable=True
            )

    common.print_start_message(args)

    _, anno, *_ = common.load_references(args, load_genome=False, \
        load_proteome=True, load_canonical_peptides=False,
        check_protein_coding=True)

    tx2gene = {}
    coding_tx = set()
    for tx_id in anno.transcripts:
        tx_model = anno.transcripts[tx_id]
        tx2gene[tx_id] = tx_model.transcript.gene_id
        if tx_model.is_protein_coding:
            coding_tx.add(tx_id)
    del anno

    if args.order_source:
        source_order = {}
        for i,val in enumerate(args.order_source.split(',')):
            if val in source_order:
                raise ValueError(
                    f"Non-unique value found from `--group-source`: {val}"
                )
            source_order[val] = i
    else:
        source_order = None

    group_map = None
    if args.group_source:
        group_map = {}
        for it in args.group_source:
            key, val = it.split(':')
            for v in val.split(','):
                group_map[v] = key

    splitter = PeptidePoolSplitter(order=source_order, group_map=group_map)

    with open(args.variant_peptides, 'rt') as handle:
        splitter.load_database(handle)

    logger(f"Variant FASTA loaded: {args.variant_peptides}.")

    if args.noncoding_peptides:
        with open(args.noncoding_peptides, 'rt') as handle:
            splitter.load_database(handle)
        logger(f"Noncoding FASTA loaded: {args.noncoding_peptides}")

    if args.alt_translation_peptides:
        with open(args.alt_translation_peptides, 'rt') as handle:
            splitter.load_database(handle)
        logger(f"Alternative Translation FASTA loaded: {args.alt_translation_peptides}")

    for file in args.gvf:
        with open(file, 'rt') as handle:
            splitter.load_gvf(handle)
        logger(f"GVF file used: {file}")

    additional_split = args.additional_split or []
    sep = SPLIT_DATABASE_KEY_SEPARATER
    additional_split = [set(x.split(sep)) for x in additional_split]

    logger(f"Using source order: {splitter.order}")
    logger(f"Using source group: {splitter.get_reversed_group_map()}")

    logger("Start splitting...")

    splitter.split(
        max_groups=args.max_source_groups,
        additional_split=additional_split,
        tx2gene=tx2gene,
        coding_tx=coding_tx
    )

    if not args.quiet:
        logger('Database split finished')

    splitter.write(args.output_prefix)

    if not args.quiet:
        logger('Split databases saved to disk.')
