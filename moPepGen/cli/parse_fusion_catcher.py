""" Module for FusionCatcher parser """
from typing import List
import pathlib
import argparse
import pickle
from moPepGen import logger, gtf, seqvar, parser, dna
from .common import add_args_reference, add_args_verbose, \
    print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_parse_fusion_catcher(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseFusionCatcher """

    p = subparsers.add_parser(
        name='parseFusionCatcher',
        help='Parse FusionCatcher result for moPepGen to call variant peptides.',
        description='Parse the FusionCatcher result to TVF format of variant'
        'records for moPepGen to call variant peptides. The genome'
    )
    p.add_argument(
        '-f', '--fusion',
        type=str,
        help="Path to the FusionCatcher's output file.",
        metavar='',
        required=True
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar='',
        required=True
    )
    add_args_reference(p)
    add_args_verbose(p)
    p.set_defaults(func=parse_fusion_catcher)
    print_help_if_missing_args(p)

def parse_fusion_catcher(args:argparse.Namespace) -> None:
    """ Parse FusionCatcher output and save it in TVF format. """
    # unpack args
    fusion = args.fusion
    index_dir:str = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.tvf'
    verbose = args.verbose

    if verbose:
        logger('moPepGen parseFusionCatcher started.')

    if index_dir:
        with open(f'{index_dir}/genome.pickle', 'rb') as handle:
            genome:dna.DNASeqDict = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            anno:gtf.GenomicAnnotation = pickle.load(handle)

        if verbose:
            logger('Indexed genome and annotation loaded.')

    else:
        genome_fasta:str = args.genome_fasta
        annotation_gtf:str = args.annotation_gtf

        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

    variants:List[seqvar.VariantRecord] = []

    for record in parser.FusionCatcherParser.parse(fusion):
        var_records = record.convert_to_variant_records(anno, genome)
        variants.extend(var_records)

    if verbose:
        logger(f'FusionCatcher output {fusion} loaded.')

    variants.sort()

    if verbose:
        logger('Variants sorted.')

    if index_dir:
        reference_index = pathlib.Path(index_dir).absolute()
        genome_fasta = None
        annotation_gtf = None
    else:
        reference_index = None
        genome_fasta = pathlib.Path(genome_fasta).absolute()
        annotation_gtf = pathlib.Path(annotation_gtf).absolute()

    metadata = seqvar.TVFMetadata(
        parser='parseFusionCatcher',
        reference_index=reference_index,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )

    seqvar.io.write(variants, output_path, metadata)
