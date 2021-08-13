""" Module for CIRCexplorer parser """
import argparse
from typing import List, Dict
import pickle
from pathlib import Path
from moPepGen import logger, dna, gtf, circ
from moPepGen.seqvar import TVFMetadata
from moPepGen.parser import CIRCexplorerParser
from .common import add_args_reference, add_args_verbose, \
    print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_parse_circexplorer(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseCIRCexplorer """
    p = subparsers.add_parser(
        name='parseCIRCexplorer',
        help='Parse CIRCexplorer result',
        description='Parse CIRCexplorer result to a TSV format for moPepGen to'
        ' call variant peptides'
    )
    p.add_argument(
        '-i', '--input-path',
        type=Path,
        help='The input file path for CIRCexplorer result. Only the known'
        'circRNA result is supported.',
        required=True,
        metavar=''
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Output prefix',
        required=True,
        metavar=''
    )
    add_args_reference(p)
    add_args_verbose(p)
    p.set_defaults(func=parse_circexplorer)
    print_help_if_missing_args(p)

def parse_circexplorer(args:argparse.Namespace):
    """ Parse circexplorer known circRNA results. """
    input_path = args.input_path
    output_prefix = args.output_prefix
    output_path = output_prefix + '.tsv'

    index_dir:Path = args.index_dir
    genome_fasta:Path = args.genome_fasta
    annotation_gtf:Path = args.annotation_gtf
    verbose = args.verbose

    if verbose:
        logger('moPepGen parseCIRCexplorer started.')

    if index_dir:
        with open(index_dir/'genome.pickle', 'rb') as handle:
            genome:dna.DNASeqDict = pickle.load(handle)

        with open(index_dir/'annotation.pickle', 'rb') as handle:
            anno:gtf.GenomicAnnotation = pickle.load(handle)

        if verbose:
            logger('Indexed genome and annotation loaded.')

    else:
        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

    circ_records:Dict[str, List[circ.CircRNAModel]] = {}

    for record in CIRCexplorerParser.parse(input_path):
        circ_record = record.convert_to_circ_rna(anno)
        gene_id = circ_record.gene_id
        if gene_id not in circ_records:
            circ_records[gene_id] = []
        circ_records[gene_id].append(circ_record)

    records = []
    for val in circ_records.values():
        records.extend(val)

    metadata = TVFMetadata(
        parser='parseCIRCexplorer',
        reference_index=index_dir,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )
    circ.io.write(records, metadata, output_path)
