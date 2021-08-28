""" Module for calling noncoding peptide """
from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path
import pkg_resources
from moPepGen import svgraph, aa, logger
from .common import add_args_cleavage, add_args_verbose, add_args_reference, \
    print_start_message, print_help_if_missing_args, load_references


if TYPE_CHECKING:
    import argparse

# pylint: disable=W0212
def add_subparser_call_noncoding(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callNoncoding """
    p = subparsers.add_parser(
        name='callNoncoding',
        help='Call non-canonical peptides from noncoding transcripts.'
    )
    p.add_argument(
        '-t', '--min-tx-length',
        type=int,
        help='Minimal transcript length.',
        metavar='',
        default=21
    )
    p.add_argument(
        '-i', '--inclusion-biotypes',
        type=Path,
        help='Inclusion biotype list.',
        metavar='',
        default=None
    )
    p.add_argument(
        '-e', '--exclusive-biotypes',
        type=Path,
        help='Exclusion biotype list.',
        metavar='',
        default=None
    )
    p.add_argument(
        '-o', '--output-fasta',
        type=str,
        help='Filename for the output FASTA.',
        metavar='',
        required=True
    )

    add_args_reference(p)
    add_args_cleavage(p)
    add_args_verbose(p)

    p.set_defaults(func=call_noncoding_peptide)
    print_help_if_missing_args(p)

def call_noncoding_peptide(args:argparse.Namespace) -> None:
    """ Main entry poitn for calling noncoding peptide """
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    min_length:int = args.min_length
    max_length:int = args.max_length

    print_start_message(args)

    genome, anno, proteome, canonical_peptides = load_references(
        args=args, load_proteome=True
    )

    inclusion_biotypes = []
    if args.inclusion_biotypes:
        with open(args.inclusion_biotypes, 'rt') as handle:
            for line in handle:
                inclusion_biotypes.append(line.rstrip())

    exclusion_path = args.exclusion_biotypes
    if not exclusion_path:
        exclusion_path = pkg_resources.resource_filename(
            'moPepGen', 'data/gencode_hs_exclusion_list.txt'
        )

    exclusion_biotypes = []
    if exclusion_path:
        with open(exclusion_path, 'rt') as handle:
            for line in handle:
                exclusion_biotypes.append(line.rstrip())

    noncanonical_pool = aa.VariantPeptidePool()

    i = 0
    for tx_id, tx_model in anno.transcripts.items():
        if inclusion_biotypes and \
                tx_model.transcript.biotype not in inclusion_biotypes:
            continue
        if exclusion_biotypes and \
                tx_model.transcript.biotype in exclusion_biotypes:
            continue
        if tx_id in proteome:
            continue
        if tx_model.transcript_len() < args.min_length:
            continue

        chrom = tx_model.transcript.location.seqname
        tx_seq = tx_model.get_transcript_sequence(genome[chrom])

        dgraph = svgraph.TranscriptVariantGraph(
            seq=tx_seq,
            _id=tx_id,
            cds_start_nf=True
        )
        dgraph.add_null_root()
        dgraph.find_all_orfs()
        pgraph = dgraph.translate()
        pgraph.form_cleavage_graph(rule=rule, exception=exception)
        peptides = pgraph.call_variant_peptides(
            miscleavage=miscleavage,
            check_variants=False
        )

        for peptide in peptides:
            noncanonical_pool.add_peptide(peptide, canonical_peptides,
                min_mw, min_length, max_length)

        if args.verbose:
            i += 1
            if i % 1000 == 0:
                logger(f'{i} transcripts processed.')

    noncanonical_pool.write(args.output_fasta)

    if args.verbose:
        logger('Noncanonical peptide FASTA file written to disk.')
