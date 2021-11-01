""" Validate callVariant result with bruteForce """
import argparse
from contextlib import redirect_stdout
from pathlib import Path
from Bio import SeqIO
from moPepGen.cli.common import add_args_reference, print_help_if_missing_args
from moPepGen.cli import call_variant_peptide
from moPepGen import logger
from moPepGen.util import downsample_reference, brute_force


# pylint: disable=W0212
def add_subparser_validate_variant_callilng(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='validateVariantCalling',
        help='Validate the varaint peptide calling result of the graph-based'
        ' algorithm with the brute force algorithm.'
    )
    parser.add_argument(
        '-i', '--input-gvf',
        type=Path,
        help='Input GVF file.'
    )
    parser.add_argument(
        '-t', '--tx-id',
        type=str,
        help='Transcript ID'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Output dir'
    )
    add_args_reference(parser, proteome=True, index=False)
    parser.set_defaults(func=validate_variant_calling)
    print_help_if_missing_args(parser)
    return parser

def call_downsample_reference(genome:Path, anno:Path, protein:Path, tx_id:str,
        output_dir:Path):
    """ downsample reference """
    args = argparse.Namespace()
    args.genome_fasta = genome
    args.annotation_gtf = anno
    args.protein_fasta = protein
    args.tx_list = [tx_id]
    args.gene_list = None
    args.output_dir = output_dir
    args.miscleavage = 2
    args.min_mw = 500.
    downsample_reference(args)

def extract_gvf(tx_id:str, gvf_file:Path, output_path:Path):
    """ extract GVF """
    i = 0
    with open(gvf_file, 'rt') as in_handle, open(output_path, 'wt') as out_handle:
        for line in in_handle:
            if line.startswith('#'):
                out_handle.write(line)
            elif tx_id in line:
                out_handle.write(line)
                i += 0
    if i > 10:
        logger(f"{i} variants found. The brute force caller is going to take a while.")

def call_variant(gvf_file:Path, ref_dir:Path, output_fasta:Path):
    """ call variant """
    args = argparse.Namespace()
    args.index_dir = None
    args.command = 'callPeptides'
    args.input_variant = [gvf_file]
    args.genome_fasta = ref_dir/'genome.fasta'
    args.annotation_gtf = ref_dir/'annotation.gtf'
    args.proteome_fasta = ref_dir/'proteome.fasta'
    args.circ_rna_bed = None
    args.output_fasta = output_fasta
    args.verbose = True
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25
    args.max_frameshift_dist = 75
    args.max_frameshift_num = 3
    call_variant_peptide(args=args)

def call_brute_force(gvf_file:Path, ref_dir:Path, output_path):
    """ call brute force """
    args = argparse.Namespace()
    args.input_gvf = gvf_file
    args.reference_dir = ref_dir
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500

    with open(output_path, 'wt') as handle:
        with redirect_stdout(handle):
            brute_force(args)

def assert_equal(variant_fasta:Path, brute_force_txt:Path):
    """ assert equal """
    with open(variant_fasta, 'rt') as handle:
        variant_seqs = set()
        for seq in SeqIO.parse(handle, 'fasta'):
            variant_seqs.add(str(seq.seq))
    with open(brute_force_txt, 'rt') as handle:
        brute_force_seqs = set()
        for line in handle:
            brute_force_seqs.add(line.rstrip())

    if variant_seqs != brute_force_seqs:
        logger('Not equal!')
    else:
        logger('Equal!')

def validate_variant_calling(args:argparse.Namespace):
    """ main entrypoint """
    output_dir:Path = args.output_dir
    ref_dir = output_dir/'index'
    call_downsample_reference(
        genome=args.genome_fasta,
        anno=args.annotation_gtf,
        protein=args.proteome_fasta,
        tx_id=args.tx_id,
        output_dir=ref_dir
    )

    logger('Reference files downsampling completed.')

    temp_gvf = output_dir/f"{args.tx_id}.gvf"
    extract_gvf(
        tx_id=args.tx_id,
        gvf_file=args.input_gvf,
        output_path=temp_gvf
    )

    logger('Transcript variants extracted from GVF.')

    variant_fasta = output_dir/'call_variant.fasta'
    call_variant(
        gvf_file=temp_gvf,
        ref_dir=ref_dir,
        output_fasta=variant_fasta
    )

    logger('Variant peptide calling completed.')

    brute_force_txt = output_dir/'brute_force.txt'
    call_brute_force(
        gvf_file=temp_gvf,
        ref_dir=ref_dir,
        output_path=brute_force_txt
    )

    logger('Brute force completed.')

    assert_equal(variant_fasta, brute_force_txt)
