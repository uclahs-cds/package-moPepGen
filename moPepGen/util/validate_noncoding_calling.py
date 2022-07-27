""" Validate callNoncoding results with bruteForceNoncoding """
import argparse
from contextlib import redirect_stdout
from pathlib import Path
import random
from typing import IO, Set
from Bio import SeqIO
from moPepGen import logger
from moPepGen.cli import call_noncoding_peptide, common
from moPepGen.gtf import GtfIO
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen.util.validate_variant_calling import call_downsample_reference
from moPepGen.util import brute_force_noncoding


# pylint: disable=W0212
def add_subparser_validate_noncoding_calling(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='validateNoncodingCalling',
        help='Validate the noncoding peptide calling result of the graph-based'
        ' algorithm with the brute force algorithm.'
    )
    parser.add_argument(
        '--tx-id',
        type=str,
        nargs='*',
        metavar='<value>',
        help='Transcript IDs to sample from to validate. If not given, all'
        ' noncoding transcripts from the annotation GTF will be used.'
    )
    parser.add_argument(
        '--n-iter',
        type=int,
        help='Number of iterations',
        metavar='<number>',
        default=1
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Path to the output diractory.',
        metavar='<file>'
    )
    common.add_args_reference(parser, proteome=True, index=False)
    parser.set_defaults(func=validate_noncoding_calling)
    common.print_help_if_missing_args(parser)
    return parser

def get_transcript_ids(handle:IO) -> Set[str]:
    """ Get all transcript IDs. """
    _, exclusion_biotypes = common.load_inclusion_exclusion_biotypes(
        argparse.Namespace(inclusion_biotypes=None, exclusion_biotypes=None)
    )
    tx_ids:Set[str]= set()
    record:GTFSeqFeature
    for record in GtfIO.parse(handle):
        record.source = 'GENCODE'
        if record.type.lower() == 'gene':
            continue
        if record.biotype in exclusion_biotypes:
            continue
        tx_ids.add(record.transcript_id)
    return tx_ids

def call_noncoding(ref_dir:Path, output_fasta:Path):
    """ call callNoncoding """
    args = argparse.Namespace()
    args.command = 'callNoncoding'
    args.genome_fasta = ref_dir/'genome.fasta'
    args.annotation_gtf = ref_dir/'annotation.gtf'
    args.proteome_fasta = ref_dir/'proteome.fasta'
    args.index_dir = None
    args.min_tx_length = 21
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25
    args.output_orf = None
    args.output_path = output_fasta
    args.min_tx_length = 21
    args.inclusion_biotypes = None
    args.exclusion_biotypes = None
    args.quiet = True
    call_noncoding_peptide(args)

def call_brute_force_noncoding(tx_id:str, ref_dir:Path, output_path:Path):
    """ call bruteForceNoncoding """
    args = argparse.Namespace()
    args.tx_id = tx_id
    args.reference_dir = ref_dir
    args.canonical_peptides = None
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25

    with open(output_path, 'wt') as handle:
        with redirect_stdout(handle):
            brute_force_noncoding(args)

def assert_equal(noncoding_fasta:Path, brute_force_txt:Path, output_dir:Path) -> bool:
    """ Compare the results of callNoncoding and bruteForceNoncoding """
    with open(noncoding_fasta, 'rt') as handle:
        noncoding_seqs = set()
        for seq in SeqIO.parse(handle, 'fasta'):
            noncoding_seqs.add(str(seq.seq))

    with open(brute_force_txt, 'rt') as handle:
        brute_force_seqs = set()
        for line in handle:
            brute_force_seqs.add(line.rstrip())

    if noncoding_seqs != brute_force_seqs:
        noncoding_only = noncoding_seqs - brute_force_seqs
        brute_force_only = brute_force_seqs - noncoding_seqs

        if noncoding_only:
            noncoding_only_file = output_dir/'call_noncoding_only.txt'
            with open(noncoding_only_file, 'wt') as handle:
                for seq in noncoding_only:
                    handle.write(seq + '\n')

        if brute_force_only:
            brute_force_only_file = output_dir/'brute_force_only.txt'
            with open(brute_force_only_file, 'wt') as handle:
                for seq in brute_force_only:
                    handle.write(seq + '\n')

        return False
    return True

def validate_noncoding_calling(args:argparse.Namespace):
    """ main entrypoint """
    if args.tx_id:
        tx_ids = args.tx_id
    else:
        with open(args.annotation_gtf, 'rt') as handle:
            tx_ids = get_transcript_ids(handle)

    n_iter = min(args.n_iter, len(tx_ids))
    tx_ids_sampled = random.sample(tx_ids, n_iter)

    for tx_id in tx_ids_sampled:
        output_dir:Path = args.output_dir/tx_id
        output_dir.mkdir(exist_ok=True)
        ref_dir = output_dir/'index'
        ref_dir.mkdir(exist_ok=True)

        noncoding_fasta = output_dir/'call_noncoding.fasta'
        brute_force_txt = output_dir/'brute_force_noncoding.fasta'

        call_downsample_reference(
            genome=args.genome_fasta,
            anno=args.annotation_gtf,
            protein=args.proteome_fasta,
            tx_id=tx_id,
            output_dir=ref_dir
        )

        call_noncoding(
            ref_dir=ref_dir,
            output_fasta=noncoding_fasta
        )

        call_brute_force_noncoding(
            tx_id=tx_id,
            ref_dir=ref_dir,
            output_path=brute_force_txt
        )

        res = assert_equal(
            noncoding_fasta=noncoding_fasta,
            brute_force_txt=brute_force_txt,
            output_dir=output_dir
        )

        logger(f"Transcript ID: {tx_id}, {'Equal' if res else 'Not equal'}!")
