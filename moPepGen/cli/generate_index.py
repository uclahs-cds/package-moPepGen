""" `generateIndex` takes the reference genome FASTA, annotation GTF, and the
translated proteome FASTA file, and processes them so they can be read by
subsequent moPepGen commands quickly. The outputted index files also contain the
canonical peptide pool. The index files can then be used in any moPepGen
command. It is recommended to run `generateIndex` before any analysis using
moPepGen to avoid processing the reference files repeatedly and save massive
time. """
from __future__ import annotations
import argparse
from pathlib import Path
import sys
from moPepGen import dna, aa, params, get_logger
from moPepGen.index import IndexDir
from moPepGen.cli import common


# pylint: disable=W0212
def add_subparser_generate_index(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen generateIndex """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='generateIndex',
        help='Generate genome and proteome index files for moPepGen',
        description='Generate genome and proteome index files for moPepGen'
        ' parsers and peptide caller.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Ouput directory for index files.',
        metavar='<file>',
        dest='output_dir',
        required=True
    )
    p.add_argument(
        '--gtf-symlink',
        help='Create a symlink of the GTF file instead of copying it.',
        action='store_true'
    )
    p.add_argument(
        '-f', '--force',
        action='store_true',
        help='Force write data to index dir.'
    )
    common.add_args_reference(p, index=False)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)
    p.set_defaults(func=generate_index)
    common.print_help_if_missing_args(p)
    return p

def generate_index(args:argparse.Namespace):
    """ Generate index  """
    path_genome:Path = args.genome_fasta
    path_gtf:Path = args.annotation_gtf
    parth_proteome:Path = args.proteome_fasta

    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    min_length:int = int(args.min_length)
    max_length:int = int(args.max_length)
    exception = args.cleavage_exception
    invalid_protein_as_noncoding:bool = args.invalid_protein_as_noncoding

    output_dir:Path = args.output_dir

    common.print_start_message(args)
    logger = get_logger()

    output_dir.mkdir(exist_ok=True)
    index_dir = IndexDir(output_dir)
    if any(index_dir.path.iterdir()):
        if args.force:
            index_dir.wipe_canonical_peptides()
            index_dir.init_metadata()
        else:
            logger.error("Index directory already exists.")
            sys.exit(1)

    # genome fasta
    genome = dna.DNASeqDict()
    genome.dump_fasta(path_genome)
    logger.info('Genome FASTA loaded')
    index_dir.save_genome(genome)
    logger.info('Genome FASTA saved to disk.')

    # proteome
    proteome = aa.AminoAcidSeqDict()
    proteome.dump_fasta(parth_proteome, source=args.reference_source)
    logger.info('Proteome FASTA loaded.')
    index_dir.save_proteome(proteome)
    logger.info('Proteome FASTA saved to disk.')

    # annotation GTF
    anno = index_dir.save_annotation(
        file=path_gtf,
        source=args.reference_source,
        proteome=proteome,
        invalid_protein_as_noncoding=invalid_protein_as_noncoding,
        symlink=args.gtf_symlink
    )
    logger.info('Genome annotation GTF saved to disk.')

    # canoincal peptide pool
    canonical_peptides = proteome.create_unique_peptide_pool(
        anno=anno, rule=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    cleavage_params = params.CleavageParams(
        enzyme=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    logger.info('canonical peptide pool generated.')
    index_dir.save_canonical_peptides(canonical_peptides, cleavage_params)
    logger.info('canonical peptide pool saved to disk.')

    # create list of coding transcripts
    coding_tx = {tx_id for tx_id, tx_model in anno.transcripts.items()
        if tx_model.is_protein_coding}
    index_dir.save_coding_tx(coding_tx)

    # save metadata
    index_dir.metadata.source = anno.source
    index_dir.save_metadata()
    logger.info('metadata saved.')
