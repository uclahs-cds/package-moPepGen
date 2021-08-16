""" Module for the moPepGen generateIndex subcommand """
import argparse
import pickle
from moPepGen import dna, aa, gtf, logger
from .common import add_args_cleavage, add_args_reference, add_args_verbose, \
    print_help_if_missing_args, print_start_message


# pylint: disable=W0212
def add_subparser_generate_index(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen generateIndex """
    p = subparsers.add_parser(
        name='generateIndex',
        help='Generate genome and proteome index files for moPepGen',
        description='Generate genome and proteome index files for moPepGen'
        'parsers and peptide caller.'
    )
    p.add_argument(
        '-o', '--output-dir',
        type=str,
        help='Ouput directory for index files.',
        metavar='',
        dest='output_dir',
        required=True
    )
    add_args_reference(p)
    add_args_cleavage(p)
    add_args_verbose(p)
    p.set_defaults(func=generate_index)
    print_help_if_missing_args(p)


def generate_index(args:argparse.Namespace):
    """ Generate  """
    path_genome:str = args.genome_fasta
    path_gtf:str = args.annotation_gtf
    parth_proteome:str = args.proteome_fasta

    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    min_length:int = int(args.min_length)
    max_length:int = int(args.max_length)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    verbose:bool = args.verbose

    output_dir:str = args.output_dir
    output_genome = f"{output_dir}/genome.pickle"
    output_proteome = f"{output_dir}/proteome.pickle"
    output_anno = f"{output_dir}/annotation.pickle"
    output_peptides = f"{output_dir}/canonical_peptides.pickle"

    print_start_message(args)

    genome = dna.DNASeqDict()
    genome.dump_fasta(path_genome)
    if verbose:
        logger('Genome FASTA loaded')
    with open(output_genome, 'wb') as handle:
        pickle.dump(genome, handle)
    if verbose:
        logger('Genome FASTA saved to disk.')
    del genome

    anno = gtf.GenomicAnnotation()
    anno.dump_gtf(path_gtf)
    if verbose:
        logger('Genome annotation GTF loaded.')
    with open(output_anno, 'wb') as handle:
        pickle.dump(anno, handle)
    if verbose:
        logger('Genome annotation GTF saved to disk.')

    proteome = aa.AminoAcidSeqDict()
    proteome.dump_fasta(parth_proteome)
    if verbose:
        logger('Proteome FASTA loaded.')
    with open(output_proteome, 'wb') as handle:
        pickle.dump(proteome, handle)
    if verbose:
        logger('Proteome FASTA saved to disk.')

    canonical_peptides = proteome.create_unique_peptide_pool(
        anno=anno, rule=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    if verbose:
        logger('canonical peptide pool generated.')
    with open(output_peptides, 'wb') as handle:
        pickle.dump(canonical_peptides, handle)
    if verbose:
        logger('canonical peptide pool saved to disk.')
