""" Common functions for cli """
import argparse
import sys
from typing import Tuple, Set
from pathlib import Path
import pickle
from moPepGen import aa, dna, gtf, logger


def print_help_if_missing_args(parser:argparse.ArgumentParser):
    """ If no args are provided, print help and exit """
    if len(sys.argv) == 2 and sys.argv[1] == parser.prog.split(' ')[1]:
        parser.print_help(sys.stderr)
        sys.exit(1)

def add_args_reference(parser:argparse.ArgumentParser):
    """ Add args for reference files including reference genome assembly FASTA,
    annotation GTF, proteome FASTA, and indexed files. """
    group = parser.add_argument_group('Reference Files')
    group.add_argument(
        '-g', '--genome-fasta',
        type=str,
        help='Path to the genome assembly FASTA file.',
        metavar=''
    )
    group.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.',
        metavar=''
    )
    group.add_argument(
        '-p', '--proteome-fasta',
        type=str,
        help='Path to the translated protein sequence FASTA file. Must come'
        'from ENSEMBL/GENCODE with the same version of the genome FASTA.',
        metavar=''
    )
    group.add_argument(
        '--index-dir',
        type=str,
        help='Path to the directory of index files generated by moPepGen'
        'generateIndex. If given, --geome-fasta, --proteome-fasta and'
        '--anntotation-gtf will be ignored.',
        metavar='',
        nargs='?',
        default=None
    )

def add_args_cleavage(parser:argparse.ArgumentParser):
    """ Add args for cleavage """
    group = parser.add_argument_group('Cleavage Parameters')
    group.add_argument(
        '-c', '--cleavage-rule',
        type=str,
        help='Cleavage rule. Defaults to trypsin. Defaults to trypsin',
        default='trypsin',
        metavar=''
    )
    group.add_argument(
        '-m', '--miscleavage',
        type=int,
        help='Number of cleavages to allow. Defaults to 2',
        default=2,
        metavar=''
    )
    group.add_argument(
        '-w', '--min-mw',
        type=float,
        help='The minimal molecular weight of the non-canonical peptides.'
        'Defaults to 500',
        default=500.,
        metavar=''
    )
    group.add_argument(
        '-l', '--min-length',
        type=int,
        help='The minimal length of non-canonical peptides, inclusive.'
        'Defaults to 7',
        default=7,
        metavar=''
    )
    group.add_argument(
        '-x', '--max-length',
        type=int,
        help='The maximum length of non-canonical peptides, inclusive.'
        'Defaults to 25',
        default=25,
        metavar=''
    )

def add_args_verbose(parser:argparse.ArgumentParser):
    """ Add verbose """
    parser.add_argument(
        '-v', '--verbose',
        type=bool,
        help='Verbose',
        metavar='',
        default=True
    )


def load_references(args:argparse.Namespace, load_canonical_peptides:bool=True,
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation, Set[str]]:
    """"""
    index_dir:Path = args.index_dir
    genome_fasta:Path = args.genome_fasta
    annotation_gtf:Path = args.annotation_gtf
    proteome_fasta:Path = args.proteome_fasta
    verbose:bool = args.verbose

    if verbose:
        logger('moPepGen callPeptide started.')

    canonical_peptides = None
    if index_dir:
        with open(f'{index_dir}/genome.pickle', 'rb') as handle:
            genome = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            annotation = pickle.load(handle)

        if load_canonical_peptides:
            with open(f"{index_dir}/canonical_peptides.pickle", 'rb') as handle:
                canonical_peptides = pickle.load(handle)
    else:
        annotation = gtf.GenomicAnnotation()
        annotation.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

        proteome = aa.AminoAcidSeqDict()
        proteome.dump_fasta(proteome_fasta)
        if verbose:
            logger('Proteome FASTA loaded.')

        if load_canonical_peptides:
            rule:str = args.cleavage_rule
            miscleavage:int = int(args.miscleavage)
            min_mw:float = float(args.min_mw)
            exception = 'trypsin_exception' if rule == 'trypsin' else None
            min_length:int = args.min_length
            max_length:int = args.max_length
            canonical_peptides = proteome.create_unique_peptide_pool(
                anno=annotation, rule=rule, exception=exception,
                miscleavage=miscleavage, min_mw=min_mw, min_length=min_length,
                max_length=max_length
            )
            if verbose:
                logger('canonical peptide pool generated.')

    return genome, annotation, canonical_peptides
