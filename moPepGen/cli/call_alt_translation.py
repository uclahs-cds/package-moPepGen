""" `callAltTranslation` calls peptide sequences from coding transcripts that
harbor any alternative translation event. """
from __future__ import annotations
import argparse
from pathlib import Path
from moPepGen import params, aa, get_logger
from moPepGen.cli import common
from moPepGen.pipeline.call_alt_translation_worker import call_alt_translation_for_transcript

OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_call_alt_translation(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callAltTranslation """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callAltTranslation',
        help='Call peptides with alternative translation from coding transcripts.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help='Output path to the alternative translation peptide FASTA.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--w2f-reassignment',
        action='store_true',
        help='Include peptides with W > F (Tryptophan to Phenylalanine) '
        'reassignment.'
    )
    p.add_argument(
        '--selenocysteine-termination',
        action='store_true',
        help='Include peptides of selenoprotiens that the UGA is treated as '
        'termination instead of Sec.'
    )

    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)

    p.set_defaults(func=call_alt_translation)
    common.print_help_if_missing_args(p)
    return p


def call_alt_translation(args: argparse.Namespace) -> None:
    """Main entrypoint for calling alternative translation peptides."""
    logger = get_logger()

    common.validate_file_format(
        args.output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    cleavage_params = params.CleavageParams(
        enzyme=args.cleavage_rule,
        exception=args.cleavage_exception,
        miscleavage=int(args.miscleavage),
        min_mw=float(args.min_mw),
        min_length=args.min_length,
        max_length=args.max_length
    )

    if not (args.selenocysteine_termination or args.w2f_reassignment):
        raise ValueError(
            'At least one of --selenocysteine-termination and --w2f-reassignment'
            ' must be given.'
        )

    common.print_start_message(args)

    # Load references in CLI layer
    ref_data = common.load_references(
        args=args, load_proteome=True, cleavage_params=cleavage_params,
        load_codon_tables=True
    )

    peptide_pool = aa.VariantPeptidePool()

    # Process each coding transcript
    for tx_id in ref_data.anno.transcripts:
        tx_model = ref_data.anno.transcripts[tx_id]
        if not tx_model.is_protein_coding:
            continue

        codon_table = ref_data.codon_tables[tx_model.transcript.chrom]

        try:
            peptides = call_alt_translation_for_transcript(
                tx_id=tx_id,
                tx_model=tx_model,
                genome=ref_data.genome,
                anno=ref_data.anno,
                codon_table=codon_table,
                cleavage_params=cleavage_params,
                w2f_reassignment=args.w2f_reassignment,
                sec_truncation=args.selenocysteine_termination
            )
        except:
            logger.error('Exception raised from %s', tx_id)
            raise

        for peptide in peptides:
            peptide_pool.add_peptide(
                peptide=peptide,
                canonical_peptides=ref_data.canonical_peptides,
                cleavage_params=cleavage_params
            )

    peptide_pool.write(args.output_path)

    logger.info('Alternative translation peptide FASTA file written to disk.')
