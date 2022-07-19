""" Validate callVariant result with bruteForce """
import argparse
from contextlib import redirect_stdout
from pathlib import Path
from typing import List
import sys
from Bio import SeqIO
from moPepGen.cli import call_variant_peptide, common
from moPepGen import logger
from moPepGen.util import downsample_reference, brute_force


INPUT_FILE_FORMATS = ['.gvf']

# pylint: disable=W0212
def add_subparser_validate_variant_callilng(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='validateVariantCalling',
        help='Validate the varaint peptide calling result of the graph-based'
        ' algorithm with the brute force algorithm.'
    )
    common.add_args_input_path(
        parser, INPUT_FILE_FORMATS,
        plural=True,
        message='File path to the GVF files generated by moPepGen parsers.'
    )
    parser.add_argument(
        '--tx-id',
        type=str,
        help='Transcript ID',
        metavar='<value>'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Path to the output diractory.',
        metavar='<file>'
    )
    parser.add_argument(
        '-f', '--force',
        action='store_true',
        help='If not set, the program stops when there are more than 10'
        ' variants. When this flag is set, the program runs anyway. Noted that '
        ' the runtime is going to increase quickly after 10 variants.',
        default=False
    )
    parser.add_argument(
        '--variant-ids',
        type=str,
        help='List of variant labels.',
        nargs='*',
        metavar='<values>'
    )
    common.add_args_reference(parser, proteome=True, index=False)
    parser.set_defaults(func=validate_variant_calling)
    common.print_help_if_missing_args(parser)
    return parser

def call_downsample_reference(genome:Path, anno:Path, protein:Path, tx_id:str,
        output_dir:Path):
    """ downsample reference """
    args = argparse.Namespace()
    args.genome_fasta = genome
    args.annotation_gtf = anno
    args.proteome_fasta = protein
    args.tx_list = [tx_id]
    args.gene_list = None
    args.output_dir = output_dir
    args.miscleavage = 2
    args.min_mw = 500.
    args.translate_noncoding = 'true'
    downsample_reference(args)

def extract_gvf(tx_id:str, gvf_files:List[Path], output_dir:Path) -> List[Path]:
    """ extract GVF """
    i = 0
    temp_gvfs = []
    for j, gvf_file in enumerate(gvf_files):
        output_path = output_dir/f"{tx_id}_{j}_{gvf_file.name}"
        temp_gvfs.append(output_path)
        with open(gvf_file, 'rt') as in_handle, open(output_path, 'wt') as out_handle:
            for line in in_handle:
                if line.startswith('#'):
                    out_handle.write(line)
                elif tx_id in line:
                    out_handle.write(line)
                    i += 0
    if i > 10:
        logger(f"{i} variants found. The brute force caller is going to take a while.")
    return temp_gvfs

def call_variant(gvf_files:Path, ref_dir:Path, output_fasta:Path):
    """ call variant """
    args = argparse.Namespace()
    args.index_dir = None
    args.command = 'callPeptides'
    args.input_path = gvf_files
    args.genome_fasta = ref_dir/'genome.fasta'
    args.annotation_gtf = ref_dir/'annotation.gtf'
    args.proteome_fasta = ref_dir/'proteome.fasta'
    args.invalid_protein_as_noncoding = False
    args.output_path = output_fasta
    args.quiet = False
    args.verbose_level = 1
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25
    args.max_variants_per_node = 7
    args.additional_variants_per_misc = 2
    args.min_nodes_to_collapse = 30
    args.naa_to_collapse = 5
    args.noncanonical_transcripts = False
    args.threads = 1
    call_variant_peptide(args=args)

def call_brute_force(gvf_files:Path, ref_dir:Path, output_path:str, force:bool,
        variant_ids:List[str]):
    """ call brute force """
    args = argparse.Namespace()
    args.input_gvf = gvf_files
    args.reference_dir = ref_dir
    args.force = force
    args.variant_ids = variant_ids
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25

    with open(output_path, 'wt') as handle:
        with redirect_stdout(handle):
            brute_force(args)

def assert_equal(variant_fasta:Path, brute_force_txt:Path, output_dir:Path):
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
        variant_only = variant_seqs - brute_force_seqs
        brute_force_only = brute_force_seqs - variant_seqs
        if variant_only:
            variant_only_file = output_dir/'call_variant_only.txt'
            with open(variant_only_file, 'wt') as handle:
                for seq in variant_only:
                    handle.write(seq + '\n')
        if brute_force_only:
            brute_force_only_file = output_dir/'brute_force_only.txt'
            with open(brute_force_only_file, 'wt') as handle:
                for seq in brute_force_only:
                    handle.write(seq + '\n')
        sys.exit(1)
    logger('Equal!')

def validate_variant_calling(args:argparse.Namespace):
    """ main entrypoint """
    for file in args.input_path:
        common.validate_file_format(file, INPUT_FILE_FORMATS, True)
    output_dir:Path = args.output_dir
    output_dir.mkdir(exist_ok=True)
    ref_dir = output_dir/'index'
    ref_dir.mkdir(exist_ok=True)
    call_downsample_reference(
        genome=args.genome_fasta,
        anno=args.annotation_gtf,
        protein=args.proteome_fasta,
        tx_id=args.tx_id,
        output_dir=ref_dir
    )

    logger('Reference files downsampling completed.')

    temp_gvfs = extract_gvf(
        tx_id=args.tx_id,
        gvf_files=args.input_path,
        output_dir=output_dir
    )

    logger('Transcript variants extracted from GVF.')

    variant_fasta = output_dir/'call_variant.fasta'
    call_variant(
        gvf_files=temp_gvfs,
        ref_dir=ref_dir,
        output_fasta=variant_fasta
    )

    logger('Variant peptide calling completed.')

    brute_force_txt = output_dir/'brute_force.txt'
    call_brute_force(
        gvf_files=temp_gvfs,
        ref_dir=ref_dir,
        output_path=brute_force_txt,
        force=args.force,
        variant_ids=args.variant_ids
    )

    logger('Brute force completed.')

    assert_equal(variant_fasta, brute_force_txt, output_dir)
