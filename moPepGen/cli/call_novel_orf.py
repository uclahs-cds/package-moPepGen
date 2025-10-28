""" `callNovelORF` calls noncanonical peptide sequences from novel ORFs.
It finds all start codons of any novel ORF gene. """
from __future__ import annotations
import argparse
from pathlib import Path
from moPepGen import params, aa, get_logger
from moPepGen.err import ReferenceSeqnameNotFoundError
from moPepGen.cli import common
from moPepGen.pipeline.call_novel_orf_worker import (
    call_novel_orf_for_transcript,
    write_orf_fasta
)

OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_call_novel_orf(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callNovelORF """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callNovelORF',
        help='Call non-canonical peptides from novel ORFs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help='Output path to the novel ORF peptide FASTA.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--output-orf',
        type=Path,
        help='Output path to the FASTA file with novel ORF sequences.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=False,
        default=None
    )
    p.add_argument(
        '--min-tx-length',
        type=int,
        help='Minimal transcript length.',
        metavar='<number>',
        default=21
    )
    p.add_argument(
        '--orf-assignment',
        type=str,
        help='Defines how ORF assignment should be done. The last ORF upstream'
        ' to the peptide is used for `max` and the first (most upstream) one is'
        ' used for `min`',
        choices=['max', 'min'],
        default='max',
        metavar='<choice>'
    )
    p.add_argument(
        '--coding-novel-orf',
        action='store_true',
        help='Include coding transcripts to find alternative ORFs.'
    )
    p.add_argument(
        '--w2f-reassignment',
        action='store_true',
        help='Include peptides with W > F (Tryptophan to Phenylalanine) '
        'reassignment.'
    )
    p.add_argument(
        '--inclusion-biotypes',
        type=Path,
        help='Inclusion biotype list.',
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--exclusion-biotypes',
        type=Path,
        help='Exclusion biotype list.',
        metavar='<file>',
        default=None
    )

    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)

    p.set_defaults(func=call_novel_orf_peptide)
    common.print_help_if_missing_args(p)
    return p

def call_novel_orf_peptide(args: argparse.Namespace) -> None:
    """Main entrypoint for calling novel ORF peptides."""
    logger = get_logger()

    common.validate_file_format(
        args.output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )
    if args.output_orf:
        common.validate_file_format(
            args.output_orf, OUTPUT_FILE_FORMATS, check_writable=True
        )

    cleavage_params = params.CleavageParams(
        enzyme=args.cleavage_rule,
        exception=args.cleavage_exception,
        miscleavage=int(args.miscleavage),
        min_mw=float(args.min_mw),
        min_length=args.min_length,
        max_length=args.max_length
    )

    common.print_start_message(args)

    # Load references in CLI layer
    ref_data = common.load_references(
        args=args, load_proteome=True, cleavage_params=cleavage_params,
        load_codon_tables=True
    )

    inclusion_biotypes, exclusion_biotypes = common.load_inclusion_exclusion_biotypes(args)

    novel_orf_peptide_pool = aa.VariantPeptidePool()
    orf_pool = []

    i = 0
    for tx_id in ref_data.anno.transcripts:
        tx_model = ref_data.anno.transcripts[tx_id]
        codon_table = ref_data.codon_tables[tx_model.transcript.chrom]

        # Filter transcripts
        if tx_model.is_protein_coding:
            if not args.coding_novel_orf:
                continue
        else:
            if inclusion_biotypes and \
                    tx_model.transcript.biotype not in inclusion_biotypes:
                continue
            if exclusion_biotypes and \
                    tx_model.transcript.biotype in exclusion_biotypes:
                continue
            if tx_id in ref_data.proteome:
                continue
            if tx_model.transcript_len() < args.min_tx_length:
                continue

        try:
            peptides, orfs = call_novel_orf_for_transcript(
                tx_id=tx_id,
                tx_model=tx_model,
                genome=ref_data.genome,
                canonical_peptides=ref_data.canonical_peptides,
                codon_table=codon_table,
                cleavage_params=cleavage_params,
                orf_assignment=args.orf_assignment,
                w2f_reassignment=args.w2f_reassignment
            )

            if not orfs:
                continue

            orf_pool.extend(orfs)

            for peptide in peptides:
                novel_orf_peptide_pool.add_peptide(
                    peptide=peptide,
                    canonical_peptides=ref_data.canonical_peptides,
                    cleavage_params=cleavage_params
                )
        except ReferenceSeqnameNotFoundError as e:
            if not ReferenceSeqnameNotFoundError.raised:
                logger.warning('%s: Make sure your GTF and FASTA files match.', e.args[0])
                ReferenceSeqnameNotFoundError.mute()
        except:
            logger.error('Exception raised from %s', tx_id)
            raise

        i += 1
        if i % 5000 == 0:
            logger.info('%i transcripts processed.', i)

    novel_orf_peptide_pool.write(args.output_path)
    if args.output_orf:
        with open(args.output_orf, 'w') as handle:
            write_orf_fasta(orf_pool, handle)

    logger.info('Noncanonical peptide FASTA file written to disk.')

