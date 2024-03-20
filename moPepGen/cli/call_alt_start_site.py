""" `callAltStart` calls peptide sequences from coding transcript with alternative
start site. """
from __future__ import annotations
import argparse
from typing import TYPE_CHECKING
from pathlib import Path
from moPepGen import params, svgraph, aa, get_logger, VARIANT_PEPTIDE_SOURCE_DELIMITER
from moPepGen.err import ReferenceSeqnameNotFoundError
from moPepGen.cli import common
from moPepGen.cli.call_noncoding_peptide import get_orf_sequences, write_orf


if TYPE_CHECKING:
    from typing import Set, Tuple, List, Dict
    from Bio.Seq import Seq
    from moPepGen.dna import DNASeqDict
    from moPepGen.gtf import TranscriptAnnotationModel

OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_call_alt_start(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callAltStart """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callAltStart',
        help='Call peptides from coding transcript with alternative start site.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help='Output path to alternative start peptide FASTA.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--output-orf',
        type=Path,
        help='Output path to the FASTA file with alternative start site ORF sequences.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=False,
        default=None
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
        '--w2f-reassignment',
        action='store_true',
        help='Include peptides with W > F (Tryptophan to Phenylalanine) '
        'reassignment.'
    )

    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)

    p.set_defaults(func=call_alt_start)
    common.print_help_if_missing_args(p)
    return p

def call_alt_start(args:argparse.Namespace) -> None:
    """ Main entryponit for calling alternative start peptides. """
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

    genome, anno, proteome, canonical_peptides = common.load_references(
        args=args, load_proteome=True, cleavage_params=cleavage_params
    )

    alt_start_pool = aa.VariantPeptidePool()
    orf_pool = []

    i = 0
    for tx_id in anno.transcripts:
        tx_model = anno.transcripts[tx_id]
        if tx_id not in proteome:
            continue

        try:
            peptides, orfs = call_alt_start_peptide_main(
                tx_id=tx_id, tx_model=tx_model, genome=genome,
                canonical_peptides=canonical_peptides,
                cleavage_params=cleavage_params,
                orf_assignment=args.orf_assignment,
                w2f_reassignment=args.w2f_reassignment
            )

            if not orfs:
                continue

            orf_pool.extend(orfs)

            for peptide in peptides:
                alt_start_pool.add_peptide(peptide, canonical_peptides,
                    cleavage_params)
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

    alt_start_pool.write(args.output_path)
    if args.output_orf:
        with open(args.output_orf, 'w') as handle:
            write_orf(orf_pool, handle)

    logger.info('Noncanonical peptide FASTA file written to disk.')

def call_alt_start_peptide_main(tx_id:str, tx_model:TranscriptAnnotationModel,
        genome:DNASeqDict, canonical_peptides:Set[str],
        cleavage_params:params.CleavageParams, orf_assignment:str,
        w2f_reassignment:bool
        ) -> Tuple[Set[aa.AminoAcidSeqRecord],List[aa.AminoAcidSeqRecord]]:
    """ call peptides from alternative start site """
    chrom = tx_model.transcript.location.seqname
    try:
        contig_seq = genome[chrom]
    except KeyError as e:
        raise ReferenceSeqnameNotFoundError(chrom) from e

    tx_seq = tx_model.get_transcript_sequence(contig_seq)

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=True,
        has_known_orf=False,
        cleavage_params=cleavage_params,
        gene_id=tx_model.gene_id,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
    )
    dgraph.init_three_frames()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph()
    peptide_anno = pgraph.call_variant_peptides(
        check_variants=False,
        check_orf=True,
        denylist=canonical_peptides,
        orf_assignment=orf_assignment,
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
    orfs = get_orf_sequences(
        pgraph=pgraph,
        tx_id=tx_id,
        gene_id=tx_model.gene_id,
        tx_seq=tx_seq,
        exclude_canonical_orf=True
    )
    return peptides, orfs
