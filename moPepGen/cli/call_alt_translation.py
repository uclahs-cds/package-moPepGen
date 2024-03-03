""" `callAltTranslation` calls peptide sequences from coding transcripts that
harbor any alternative translation event. """
from __future__ import annotations
import argparse
from typing import TYPE_CHECKING, Dict, Set
from pathlib import Path
from moPepGen import params, svgraph, aa, get_logger, VARIANT_PEPTIDE_SOURCE_DELIMITER
from moPepGen.cli import common


if TYPE_CHECKING:
    from Bio.Seq import Seq
    from moPepGen.gtf import TranscriptAnnotationModel, GenomicAnnotation
    from moPepGen.dna import DNASeqDict

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


def call_alt_translation(args:argparse.Namespace) -> None:
    """ Main entrypoint for calling alternative translation peptides """
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

    genome, anno, _, canonical_peptides = common.load_references(
        args=args, load_proteome=True, cleavage_params=cleavage_params
    )

    peptide_pool = aa.VariantPeptidePool()

    for tx_id in anno.transcripts:
        tx_model = anno.transcripts[tx_id]
        if not tx_model.is_protein_coding:
            continue

        try:
            peptides = call_alt_translation_main(
                tx_id=tx_id, tx_model=tx_model,
                genome=genome, anno=anno,
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
                canonical_peptides=canonical_peptides,
                cleavage_params=cleavage_params
            )

    peptide_pool.write(args.output_path)

    logger.info('Alternative translation peptide FASTA file written to disk.')

def call_alt_translation_main(tx_id:str, tx_model:TranscriptAnnotationModel,
        genome:DNASeqDict, anno:GenomicAnnotation,
        cleavage_params:params.CleavageParams,
        w2f_reassignment:bool, sec_truncation:bool):
    """ wrapper of graph operations to call peptides """
    chrom = tx_model.transcript.chrom
    tx_seq = tx_model.get_transcript_sequence(genome[chrom])

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=True,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        cleavage_params=cleavage_params,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
    )
    dgraph.gather_sect_variants(anno)
    dgraph.init_three_frames()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph()
    peptide_anno = pgraph.call_variant_peptides(
        check_variants=True,
        truncate_sec=sec_truncation,
        w2f=w2f_reassignment,
        check_external_variants=False
    )
    peptide_map:Dict[Seq, Set[str]] = {}
    for seq, annotated_labels in peptide_anno.items():
        for label in annotated_labels:
            if seq in peptide_map:
                peptide_map[seq].add(label.label)
            else:
                peptide_map[seq] = {label.label}
    peptides = set()
    for seq, labels in peptide_map.items():
        label = VARIANT_PEPTIDE_SOURCE_DELIMITER.join(labels)
        peptides.add(
            aa.AminoAcidSeqRecord(
                seq=seq,
                description=label,
                name=label
            )
        )
    return peptides
