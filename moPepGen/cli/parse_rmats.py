""" Module for rMATS parser """
from typing import Dict, Set
import argparse
import pickle
from pathlib import Path
from moPepGen import logger, dna, gtf, seqvar
from moPepGen.parser import RMATSParser
from .common import add_args_reference, add_args_verbose, \
    print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_parse_rmats(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseRMATs """
    ## parser_rmats
    p = subparsers.add_parser(
        name='parseRMATS',
        help='Parse rMATS result for moPepGen to call variant peptides.',
        description='Parse the rMATS result to TVF format of variant'
        'records for moPepGen to call variant peptides.'
    )

    p.add_argument(
        '--skipped-exon',
        type=Path,
        help="Skipped exon junction count txt file.",
        metavar='',
        default=None,
        dest='skipped_exon'
    )
    p.add_argument(
        '--alternative-5-splicing',
        type=Path,
        help="Alternative 5' splicing junction count txt file.",
        metavar='',
        default=None,
        dest='alternative_5_splicing'
    )
    p.add_argument(
        '--alternative-3-splicing',
        type=Path,
        help="Alternative 3' splicing junction count txt file.",
        metavar='',
        default=None,
        dest='alternative_3_splicing'
    )
    p.add_argument(
        '--mutually-exclusive-exons',
        type=Path,
        help="Mutually exclusive junction count txt file.",
        metavar='',
        default=None,
        dest='mutually_exclusive_exons'
    )
    p.add_argument(
        '--retained-intron',
        type=Path,
        help="Retained intron junction count txt file.",
        metavar='',
        default=None,
        dest='retained_intron'
    )

    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar=''
    )

    add_args_reference(p)
    add_args_verbose(p)
    p.set_defaults(func=parse_rmats)
    print_help_if_missing_args(p)

def parse_rmats(args:argparse.Namespace) -> None:
    """ Parse rMATS results into TSV """
    skipped_exon = args.skipped_exon
    alternative_5 = args.alternative_5_splicing
    alternative_3 = args.alternative_3_splicing
    mutually_exclusive = args.mutually_exclusive_exons
    retained_intron = args.retained_intron
    index_dir:Path = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.tvf'
    verbose = args.verbose

    if verbose:
        logger('moPepGen parserRMATS started.')

    if index_dir:
        with open(index_dir/'genome.pickle', 'rb') as handle:
            genome:dna.DNASeqDict = pickle.load(handle)

        with open(index_dir/'annotation.pickle', 'rb') as handle:
            anno:gtf.GenomicAnnotation = pickle.load(handle)

        if verbose:
            logger('Indexed genome and annotation loaded.')

    else:
        genome_fasta:Path = args.genome_fasta
        annotation_gtf:Path = args.annotation_gtf

        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

    variants:Dict[str,Set[seqvar.VariantRecord]] = {}
    rmats_outputs = [
        ('SE', skipped_exon), ('A5SS', alternative_5), ('A3SS', alternative_3),
        ('MXE', mutually_exclusive), ('RI', retained_intron)
    ]
    for event_type, path in rmats_outputs:
        if path:
            for record in RMATSParser.parse(path, event_type):
                var_records = record.convert_to_variant_records(anno, genome)
                for var_record in var_records:
                    transcript_id = var_record.location.seqname
                    if transcript_id not in variants:
                        variants[transcript_id] = set()
                    variants[transcript_id].add(var_record)

    variants_sorted = []
    for val in variants.values():
        val = list(val)
        val.sort()
        variants_sorted.extend(val)

    if index_dir:
        reference_index = Path(index_dir).absolute()
        genome_fasta = None
        annotation_gtf = None
    else:
        reference_index = None
        genome_fasta = Path(genome_fasta).absolute()
        annotation_gtf = Path(annotation_gtf).absolute()

    metadata = seqvar.TVFMetadata(
        parser='parseRMATS',
        reference_index=reference_index,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )
    seqvar.io.write(variants_sorted, output_path, metadata)


if __name__ == '__main__':
    data = argparse.Namespace()
    data.skipped_exon = Path('test/files/alternative_splicing/rmats_se.txt')
    data.alternative_5_splicing = Path('test/files/alternative_splicing/rmats_a5ss.txt')
    data.alternative_3_splicing = Path('test/files/alternative_splicing/rmats_a3ss.txt')
    data.mutually_exclusive_exons = Path('test/files/alternative_splicing/rmats_mxe.txt')
    data.retained_intron = Path('test/files/alternative_splicing/rmats_ri.txt')
    data.index_dir = Path('test/files/index')
    data.output_prefix = 'test/files/alternative_splicing/rmats'
    data.verbose = True
    parse_rmats(data)
