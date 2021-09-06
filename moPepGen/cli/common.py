""" Common functions for cli """
import argparse
import sys
from typing import Tuple, Set
from pathlib import Path
import pickle
from moPepGen import aa, dna, gtf, logger, seqvar


def print_help_if_missing_args(parser:argparse.ArgumentParser):
    """ If no args are provided, print help and exit """
    if len(sys.argv) == 2 and sys.argv[1] == parser.prog.split(' ')[-1]:
        parser.print_help(sys.stderr)
        sys.exit(1)

def print_start_message(args:argparse.Namespace):
    """ Print the program start message """
    if args.verbose:
        logger(f'moPepGen {args.command} started')

def add_args_reference(parser:argparse.ArgumentParser, genome:bool=True,
        proteome:bool=True, index:bool=True):
    """ Add args for reference files including reference genome assembly FASTA,
    annotation GTF, proteome FASTA, and indexed files. """
    group = parser.add_argument_group('Reference Files')
    if genome:
        group.add_argument(
            '-g', '--genome-fasta',
            type=Path,
            help='Path to the genome assembly FASTA file.',
            metavar='',
            default=None
        )
    group.add_argument(
        '-a', '--annotation-gtf',
        type=Path,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.',
        metavar='',
        default=None
    )
    if proteome:
        group.add_argument(
            '-p', '--proteome-fasta',
            type=Path,
            help='Path to the translated protein sequence FASTA file. Must come'
            'from ENSEMBL/GENCODE with the same version of the genome FASTA.',
            metavar='',
            default=None
        )
    if index:
        group.add_argument(
            '--index-dir',
            type=Path,
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

def add_args_source(parser:argparse.ArgumentParser):
    """ Add source """
    parser.add_argument(
        '--source',
        type=str,
        help='Variant source (e.g. gSNP, sSNV, Fusion)',
        required=True
    )

def load_references(args:argparse.Namespace, load_genome:bool=True,
        load_canonical_peptides:bool=True, load_proteome:bool=False
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation, Set[str]]:
    """ Load reference files. If index_dir is specified, data will be loaded
    from pickles, otherwise, will read from FASTA and GTF. """
    index_dir:Path = args.index_dir
    verbose:bool = args.verbose

    genome = None
    annotation = None
    canonical_peptides = None

    if index_dir:
        if load_genome:
            with open(f'{index_dir}/genome.pickle', 'rb') as handle:
                genome = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            annotation = pickle.load(handle)

        if load_proteome:
            with open(f'{index_dir}/proteome.pickle', 'rb') as handle:
                proteome = pickle.load(handle)

        if load_canonical_peptides:
            with open(f"{index_dir}/canonical_peptides.pickle", 'rb') as handle:
                canonical_peptides = pickle.load(handle)
    else:
        annotation = gtf.GenomicAnnotation()
        annotation.dump_gtf(args.annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        if load_genome:
            genome = dna.DNASeqDict()
            genome.dump_fasta(args.genome_fasta)
            if verbose:
                logger('Genome assembly FASTA loaded.')

        if load_canonical_peptides:
            proteome = aa.AminoAcidSeqDict()
            proteome.dump_fasta(args.proteome_fasta)
            if verbose:
                logger('Proteome FASTA loaded.')
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

    if not load_proteome:
        proteome = None

    return genome, annotation, proteome, canonical_peptides

def generate_metadata(args:argparse.Namespace) -> seqvar.GVFMetadata:
    """ Generate metadata """
    if args.index_dir:
        reference_index = args.index_dir.absolute()
        genome_fasta = None
        annotation_gtf = None
    else:
        reference_index = None
        genome_fasta = args.genome_fasta.absolute() if \
            hasattr(args, 'genome_fasta') else None
        annotation_gtf = args.annotation_gtf.absolute() if \
            hasattr(args, 'annotation_gtf') else None

    return seqvar.GVFMetadata(
        parser='parseVEP',
        source=args.source,
        reference_index=reference_index,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )
