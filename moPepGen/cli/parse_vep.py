""" VEP2VariantPeptides module """
from typing import Dict, List
import pathlib
import argparse
import pickle
from moPepGen.parser import VEPParser
from moPepGen import gtf, dna, seqvar, logger
from .common import add_args_reference, add_args_verbose, \
    print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_parse_vep(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseVEP """

    p = subparsers.add_parser(
        name='parseVEP',
        help='Parse VEP output for moPepGen to call variant peptides.',
        description="Parse VEP output tsv to the TVF format of variant records"
        "for moPepGen to call variant peptides. The genome assembly FASTA and"
        "annotation GTF must come from the same GENCODE/ENSEMBL version, and"
        "must the consistent with the VEP output."
    )

    p.add_argument(
        '-i', '--vep-txt',
        type=str,
        nargs='+',
        help='Path to VEP result txt file.',
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
    p.set_defaults(func=parse_vep)
    print_help_if_missing_args(p)

def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    # unpack args
    vep_files:List[str] = args.vep_txt
    index_dir:str = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.tvf'
    verbose = args.verbose
    genome_fasta = None
    annotation_gtf = None

    if verbose:
        logger('moPepGen parseVEP started.')

    # if indexed files are given, load directly
    if index_dir:
        with open(f'{index_dir}/genome.pickle', 'rb') as handle:
            genome = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            anno = pickle.load(handle)

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

    vep_records:Dict[str, List[seqvar.VariantRecord]] = {}

    for vep_file in vep_files:
        for record in VEPParser.parse(vep_file):
            transcript_id = record.feature

            if transcript_id not in vep_records.keys():
                vep_records[transcript_id] = []

            record = record.convert_to_variant_record(anno, genome)

            vep_records[transcript_id].append(record)

        if verbose:
            logger(f'VEP file {vep_file} loaded.')

    for records in vep_records.values():
        records.sort()

    if verbose:
        logger('VEP sorting done.')

    if index_dir:
        reference_index = pathlib.Path(index_dir).absolute()
        genome_fasta = None
        annotation_gtf = None
    else:
        reference_index = None
        genome_fasta = pathlib.Path(genome_fasta).absolute()
        annotation_gtf = pathlib.Path(annotation_gtf).absolute()

    metadata = seqvar.TVFMetadata(
        parser='parseVEP',
        reference_index=reference_index,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )

    all_records = []
    for records in vep_records.values():
        all_records.extend(records)

    seqvar.io.write(all_records, output_path, metadata)

    if verbose:
        logger('Variant info written to disk.')


if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.vep_txt = [
        'test/files/vep/vep_snp.txt',
        'test/files/vep/vep_indel.txt'
    ]
    test_args.index_dir = None
    test_args.genome_fasta = 'test/files/genome.fasta'
    test_args.annotation_gtf = 'test/files/annotation.gtf'
    test_args.output_prefix = 'test/files/vep/vep'
    test_args.verbose = True
    parse_vep(args=test_args)
