""" Common functions for cli """
from __future__ import annotations
import argparse
import os
import sys
from typing import Tuple, Set, List, TYPE_CHECKING
from pathlib import Path
import errno
import signal
import functools
import time
import logging
import pkg_resources
from moPepGen import aa, dna, gtf, seqvar, get_logger, constant
from moPepGen.aa.expasy_rules import EXPASY_RULES
from moPepGen.index import IndexDir


if TYPE_CHECKING:
    from moPepGen.params import CleavageParams

def print_help_if_missing_args(parser:argparse.ArgumentParser):
    """ If no args are provided, print help and exit """
    if len(sys.argv) == 2 and sys.argv[1] == parser.prog.split(' ')[-1]:
        parser.print_help(sys.stderr)
        sys.exit(1)

def print_start_message(args:argparse.Namespace):
    """ Print the program start message """
    get_logger().info('moPepGen %s started', args.command)

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
    group.add_argument(
        '--reference-source',
        type=str,
        choices=['GENCODE', 'ENSEMBL'],
        help='Source of reference genome and annotation.',
        default=None
    )
    if proteome:
        group.add_argument(
            '-p', '--proteome-fasta',
            type=Path,
            help='Path to the translated protein sequence FASTA file. Only'
            ' ENSEMBL and GENCODE are supported. Its version must be the same'
            ' as genome FASTA and annotation GTF.',
            metavar='<file>',
            default=None
        )
        group.add_argument(
            '--invalid-protein-as-noncoding',
            action='store_true',
            help='Treat any transcript that the protein sequence is invalid ('
            ' contains the * symbol) as noncoding.'
        )
    if index:
        group.add_argument(
            '--index-dir',
            type=Path,
            help='Path to the directory of index files generated by moPepGen'
            ' generateIndex. If given, --genome-fasta, --proteome-fasta and'
            ' --anntotation-gtf will be ignored.',
            metavar='<file>',
            nargs='?',
            default=None
        )

def add_args_cleavage(parser:argparse.ArgumentParser, enzyme_only:bool=False):
    """ Add args for cleavage """
    group = parser.add_argument_group('Cleavage Parameters')
    group.add_argument(
        '-c', '--cleavage-rule',
        type=str,
        help='Enzymatic cleavage rule.',
        default='trypsin',
        metavar='<value>',
        choices=list(EXPASY_RULES.keys()) + ['None']
    )
    group.add_argument(
        '--cleavage-exception',
        type=str,
        help='Enzymatic cleavage exception.',
        default='auto',
        metavar='<value>'
    )
    if enzyme_only:
        return
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
    group.add_argument(
        '--flanking-size',
        type=int,
        default=10,
        help='Flanking size for no enzymatic cleavage.',
        metavar='<number?'
    )

def add_args_decoy(parser:argparse.ArgumentParser):
    """ add decoy fasta related arguments """
    group = parser.add_argument_group('Decoy Database Parameters')
    group.add_argument(
        '--decoy-string',
        type=str,
        default='DECOY_',
        help='The decoy string that is combined with the FASTA header for decoy'
        ' sequences.',
        metavar='<value>'
    )
    group.add_argument(
        '--decoy-string-position',
        type=str,
        choices=['prefix', 'suffix'],
        help='Should the decoy string be placed at the start or end of FASTA'
        ' headers?',
        default='prefix',
        metavar='<value>'
    )

def add_args_debug_level(parser:argparse.ArgumentParser):
    """ add debug level """
    parser.add_argument(
        '--debug-level',
        type=str,
        help='Debug level.',
        default='INFO',
        metavar='<value|number>'
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Quiet',
        default=False
    )

def add_args_output_path(parser:argparse.ArgumentParser, formats:List[str]):
    """ add output file path """
    parser.add_argument(
        '-o', '--output-path',
        type=Path,
        help=f'File path to the output file. Valid formats: {formats}',
        metavar='<file>',
        required=True
    )

def add_args_input_path(parser:argparse.ArgumentParser, formats:List[str],
        plural:bool=False, message:str=None, required:bool=True):
    """ Add input file path """
    if message is None:
        message = f"File path to the input file{'s' if plural else ''}."

    message += f"{' Can take multiple files.' if plural else ''} Valid formats: " +\
        f"{formats}"

    metavar = ["<files>"] if plural else "<file>"

    if plural:
        nargs = '+' if required else '*'
    else:
        nargs = None

    parser.add_argument(
        '-i', '--input-path',
        type=Path,
        help=message,
        nargs=nargs,
        metavar=metavar,
        required=required
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
        invalid_protein_as_noncoding:bool=False, check_protein_coding:bool=False,
        cleavage_params:CleavageParams=None
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotationOnDisk, Set[str]]:
    """ Load reference files. If index_dir is specified, data will be loaded
    from pickles, otherwise, will read from FASTA and GTF. """
    logger = get_logger()
    genome = None
    anno = None
    canonical_peptides = None
    if invalid_protein_as_noncoding:
        load_proteome = True

    if args.index_dir:
        index_dir = IndexDir(args.index_dir)
        index_dir.validate_metadata()

        if load_canonical_peptides:
            canonical_peptides = index_dir.load_canonical_peptides(cleavage_params)

        if load_genome:
            genome = index_dir.load_genome()

        anno = index_dir.load_annotation()

        if load_proteome:
            proteome = index_dir.load_proteome()

        if invalid_protein_as_noncoding:
            anno.check_protein_coding(proteome, True)

        logger.info('Reference indices loaded.')
    else:
        if (check_protein_coding is True or load_canonical_peptides is True) and \
                args.proteome_fasta is None:
            raise ValueError('--proteome-fasta was not specified.')
        anno = gtf.GenomicAnnotationOnDisk()
        anno.generate_index(args.annotation_gtf, source=args.reference_source)

        logger.info('Annotation GTF loaded.')

        if load_proteome or load_canonical_peptides or check_protein_coding:
            proteome = aa.AminoAcidSeqDict()
            proteome.dump_fasta(args.proteome_fasta, source=args.reference_source)
            logger.info('Proteome FASTA loaded.')
            anno.check_protein_coding(proteome, invalid_protein_as_noncoding)

        if load_genome:
            genome = dna.DNASeqDict()
            genome.dump_fasta(args.genome_fasta)
            logger.info('Genome assembly FASTA loaded.')

        if load_canonical_peptides:
            rule:str = args.cleavage_rule
            miscleavage:int = int(args.miscleavage)
            min_mw:float = float(args.min_mw)
            exception = args.cleavage_exception
            min_length:int = args.min_length
            max_length:int = args.max_length
            canonical_peptides = proteome.create_unique_peptide_pool(
                anno=anno, rule=rule, exception=exception,
                miscleavage=miscleavage, min_mw=min_mw, min_length=min_length,
                max_length=max_length
            )
            logger.info('canonical peptide pool generated.')

    if not load_proteome:
        proteome = None

    return genome, anno, proteome, canonical_peptides

def generate_metadata(args:argparse.Namespace) -> seqvar.GVFMetadata:
    """ Generate metadata """
    if args.index_dir:
        reference_index = args.index_dir.resolve()
        genome_fasta = None
        annotation_gtf = None
    else:
        reference_index = None
        genome_fasta = args.genome_fasta.resolve() if \
            hasattr(args, 'genome_fasta') else None
        annotation_gtf = args.annotation_gtf.resolve() if \
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

def validate_file_format(file:Path, types:List[str]=None, check_readable:bool=False,
            check_writable:bool=False, is_directory:bool=False):
    """ Validate the file type """
    if types is None:
        types = []
    if not is_directory:
        suffixes = file.suffixes
        actual_suffixes = [suffixes[-1]]
        if len(suffixes) > 1:
            actual_suffixes.append(''.join(suffixes[-2:]))
        if not any(suffix in types for suffix in actual_suffixes):
            raise ValueError(
                f'The file {file} is invalid. Valid file types are {types}.'
            )
    if check_readable:
        if not file.exists():
            raise FileNotFoundError(f"File does not exist: '{file}'")
        if not os.access(file, os.R_OK):
            raise PermissionError(f"Permission denied: '{file}'")
    if check_writable:
        if file.exists():
            if not os.access(file, os.W_OK):
                raise PermissionError(f"Permission denied: '{file}'")
        else:
            if not os.access(file.parent, os.W_OK):
                raise PermissionError(f"Permission denied: '{file}'")

# pylint: disable=unused-argument
def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    """ Decorator to raise a TimeoutError if the process runs over time. """
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if 'timeout' in kwargs and kwargs['timeout']:
                seconds = kwargs['timeout']
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator

def setup_loggers(level:str):
    """ Initialize loggers for both init and run """
    # pylint: disable=protected-access
    try:
        level = int(level)
    except ValueError:
        pass
    try:
        log_level = logging._checkLevel(level)
    except ValueError:
        log_level = logging._checkLevel('INFO')

    logger = logging.getLogger(constant.PROG_NAME)
    logger.setLevel(log_level)

    formatter = logging.Formatter(
        '[ %(asctime)s ] [ %(name)s - %(levelname)5s ] %(message)s',
        "%Y-%m-%d %H:%M:%S"
    )
    logging.Formatter.converter = time.localtime

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(log_level)

    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger
