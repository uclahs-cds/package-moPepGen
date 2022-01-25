""" Common functions for cli """
from __future__ import annotations
import argparse
import sys
from typing import Tuple, Set, List
from pathlib import Path
import pickle
import pkg_resources
from moPepGen import aa, dna, gtf, logger, seqvar, err
from moPepGen.aa.expasy_rules import EXPASY_RULES
from moPepGen.dna.DNASeqDict import DNASeqDict
from moPepGen.version import MetaVersion


def print_help_if_missing_args(parser:argparse.ArgumentParser):
    """ If no args are provided, print help and exit """
    if len(sys.argv) == 2 and sys.argv[1] == parser.prog.split(' ')[-1]:
        parser.print_help(sys.stderr)
        sys.exit(1)

def print_start_message(args:argparse.Namespace):
    """ Print the program start message """
    if not args.quiet:
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
            help='Path to the genome assembly FASTA file. Only ENSEMBL and'
            ' GENCODE are supported. Its version must be the same as the'
            ' annotation GTF and proteome FASTA',
            metavar='<file>',
            default=None
        )
    group.add_argument(
        '-a', '--annotation-gtf',
        type=Path,
        help='Path to the annotation GTF file. Only ENSEMBL and GENCODE are'
        ' supported. Its version must be the same as the genome and proteome'
        ' FASTA.',
        metavar='<file>',
        default=None
    )
    if proteome:
        group.add_argument(
            '-p', '--proteome-fasta',
            type=Path,
            help='Path to the translated protein sequence FASTA file. Only'
            ' ENSEMBL and GENCODE are supported. Its verion must be the same'
            ' as genome FASTA and annotation GTF.',
            metavar='<file>',
            default=None
        )
        group.add_argument(
            '--invalid-protein-as-noncoding',
            action='store_true',
            help='Treat any transcript that the protein sequence is invalid ('
            'contains the * symbol) as noncoding.'
        )
    if index:
        group.add_argument(
            '--index-dir',
            type=Path,
            help='Path to the directory of index files generated by moPepGen'
            'generateIndex. If given, --geome-fasta, --proteome-fasta and'
            '--anntotation-gtf will be ignored.',
            metavar='<file>',
            nargs='?',
            default=None
        )

def add_args_cleavage(parser:argparse.ArgumentParser):
    """ Add args for cleavage """
    group = parser.add_argument_group('Cleavage Parameters')
    group.add_argument(
        '-c', '--cleavage-rule',
        type=str,
        help='Enzymatic cleavage rule.',
        default='trypsin',
        metavar='<value>',
        choices=list(EXPASY_RULES.keys())
    )
    group.add_argument(
        '-m', '--miscleavage',
        type=int,
        help='Number of cleavages to allow per non-canonical peptide.',
        default=2,
        metavar='<number>'
    )
    group.add_argument(
        '-w', '--min-mw',
        type=float,
        help='The minimal molecular weight of the non-canonical peptides.',
        default=500.,
        metavar='<number>'
    )
    group.add_argument(
        '-l', '--min-length',
        type=int,
        help='The minimal length of non-canonical peptides, inclusive.',
        default=7,
        metavar='<number>'
    )
    group.add_argument(
        '-x', '--max-length',
        type=int,
        help='The maximum length of non-canonical peptides, inclusive.',
        default=25,
        metavar='<number>'
    )

def add_args_quiet(parser:argparse.ArgumentParser):
    """ Add quiet """
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Quiet',
        default=False
    )

def add_args_output_prefix(parser:argparse.ArgumentParser):
    """ add output prefix """
    parser.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename. The output file will be saved as'
        ' <output_prefix>.gvf',
        metavar='<value>',
        required=True
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
        load_canonical_peptides:bool=True, load_proteome:bool=False,
        invalid_protein_as_noncoding:bool=False
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation, Set[str]]:
    """ Load reference files. If index_dir is specified, data will be loaded
    from pickles, otherwise, will read from FASTA and GTF. """
    index_dir:Path = args.index_dir
    quiet:bool = args.quiet

    genome = None
    annotation = None
    canonical_peptides = None

    version = MetaVersion()

    if index_dir:
        if load_genome:
            with open(f'{index_dir}/genome.pkl', 'rb') as handle:
                genome:DNASeqDict = pickle.load(handle)
                if not version.is_valid(genome.version):
                    raise err.IndexVersionNotMatchError(version, genome.version)

        with open(f'{index_dir}/annotation.pkl', 'rb') as handle:
            annotation:gtf.GenomicAnnotation = pickle.load(handle)
            if not version.is_valid(annotation.version):
                raise err.IndexVersionNotMatchError(version, genome.version)


        if load_proteome:
            with not open(f'{index_dir}/proteome.pkl', 'rb') as handle:
                proteome:aa.AminoAcidSeqDict = pickle.load(handle)
                if version.is_valid(proteome.version):
                    raise err.IndexVersionNotMatchError(version, genome.version)

        if load_canonical_peptides:
            with open(f"{index_dir}/canonical_peptides.pkl", 'rb') as handle:
                canonical_peptides = pickle.load(handle)
        if not quiet:
            logger('Reference indices loaded.')
    else:
        annotation = gtf.GenomicAnnotation()
        annotation.dump_gtf(args.annotation_gtf)
        if not quiet:
            logger('Annotation GTF loaded.')

        proteome = aa.AminoAcidSeqDict()
        proteome.dump_fasta(args.proteome_fasta)
        if not quiet:
            logger('Proteome FASTA loaded.')

        annotation.check_protein_coding(proteome, invalid_protein_as_noncoding)

        if load_genome:
            genome = dna.DNASeqDict()
            genome.dump_fasta(args.genome_fasta)
            if not quiet:
                logger('Genome assembly FASTA loaded.')

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
            if not quiet:
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
        parser=args.command,
        source=args.source,
        chrom='Gene ID',
        reference_index=reference_index,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )

def load_inclusion_exclusion_biotypes(args:argparse.Namespace
        ) -> Tuple[List[str], List[str]]:
    """ Load inclusion and exclusion biotypes """
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

    return inclusion_biotypes, exclusion_biotypes

def parse_range(x:str) -> Tuple[int,int]:
    """ parse range from argument """
    y = x.split(',')
    if len(y) != 2:
        raise ValueError('Range invalid')
    return tuple(int(i) for i in y)
