""" `generateIndex` takes the reference genome FASTA, annotation GTF, and the
translated proteome FASTA file, converts them to the moPepGen objects,
serializes them and saves into disk. The outputted index files also contain the
canonical peptide pool. The index files can then be used in any moPepGen
command. It is recommended to run `generateIndex` before any analysis using
moPepGen to avoid processing the reference files repeatedly and save massive
time. """
from __future__ import annotations
import argparse
from pathlib import Path
import pickle
from moPepGen import dna, aa, gtf, logger
from .common import add_args_cleavage, add_args_reference, add_args_quiet, \
    print_help_if_missing_args, print_start_message


# pylint: disable=W0212
def add_subparser_generate_index(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen generateIndex """
    p = subparsers.add_parser(
        name='generateIndex',
        help='Generate genome and proteome index files for moPepGen',
        description='Generate genome and proteome index files for moPepGen'
        'parsers and peptide caller.',
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
    add_args_reference(p, index=False)
    add_args_cleavage(p)
    add_args_quiet(p)
    p.set_defaults(func=generate_index)
    print_help_if_missing_args(p)
    return p


def generate_index(args:argparse.Namespace):
    """ Generate  """
    path_genome:Path = args.genome_fasta
    path_gtf:Path = args.annotation_gtf
    parth_proteome:Path = args.proteome_fasta

    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    min_length:int = int(args.min_length)
    max_length:int = int(args.max_length)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    quiet:bool = args.quiet

    output_dir:Path = args.output_dir
    output_genome = output_dir/"genome.pickle"
    output_proteome = output_dir/"proteome.pickle"
    output_anno = output_dir/"annotation.pickle"
    output_peptides = output_dir/"canonical_peptides.pickle"

    print_start_message(args)

    output_dir.mkdir(exist_ok=True)

    genome = dna.DNASeqDict()
    genome.dump_fasta(path_genome)
    if not quiet:
        logger('Genome FASTA loaded')
    with open(output_genome, 'wb') as handle:
        pickle.dump(genome, handle)
    if not quiet:
        logger('Genome FASTA saved to disk.')
    del genome

    anno = gtf.GenomicAnnotation()
    anno.dump_gtf(path_gtf)
    if not quiet:
        logger('Genome annotation GTF loaded.')

    proteome = aa.AminoAcidSeqDict()
    proteome.dump_fasta(parth_proteome)

    anno.check_protein_coding(proteome)

    with open(output_anno, 'wb') as handle:
        pickle.dump(anno, handle)
    if not quiet:
        logger('Genome annotation GTF saved to disk.')

    if not quiet:
        logger('Proteome FASTA loaded.')
    with open(output_proteome, 'wb') as handle:
        pickle.dump(proteome, handle)
    if not quiet:
        logger('Proteome FASTA saved to disk.')

    canonical_peptides = proteome.create_unique_peptide_pool(
        anno=anno, rule=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    if not quiet:
        logger('canonical peptide pool generated.')
    with open(output_peptides, 'wb') as handle:
        pickle.dump(canonical_peptides, handle)
    if not quiet:
        logger('canonical peptide pool saved to disk.')
