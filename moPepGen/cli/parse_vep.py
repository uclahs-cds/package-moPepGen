""" VEP2VariantPeptides module """
from typing import Dict, List
import argparse
import pickle
from moPepGen.vep import VepIO
from moPepGen import gtf, dna, seqvar, logger


def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    # unpack args
    vep_files:List[str] = args.vep_txt
    index_dir:str = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '_moPepGem.bed'
    verbose = args.verbose

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

        anno = gtf.TranscriptGTFDict()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

    vep_records:Dict[str, List[seqvar.VariantRecord]] = {}

    for vep_file in vep_files:
        for record in VepIO.parse(vep_file):
            transcript_id = record.feature

            if transcript_id not in vep_records.keys():
                vep_records[transcript_id] = []

            chrom_seqname = record.location.split(':')[0]

            transcript_seq = anno[transcript_id]\
                .get_transcript_sequence(genome[chrom_seqname])

            record = record.convert_to_variant_record(transcript_seq)

            vep_records[transcript_id].append(record)

        if verbose:
            logger(f'VEP file {vep_file} loaded.')

    for records in vep_records.values():
        records.sort()

    if verbose:
        logger('VEP sorting done.')

    with open(output_path, 'w') as handle:
        mode = 'w'
        for records in vep_records.values():
            seqvar.io.write(records, handle, mode)
            if mode == 'w':
                mode = 'a'

    if verbose:
        logger('Variant info written to disk.')


if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.vep_txt = [
        'test/files/vep_test_files/vep_snp.txt',
        'test/files/vep_test_files/vep_indel.txt'
    ]
    test_args.index_dir = None
    test_args.genome_fasta = 'test/files/vep_test_files/genome.fasta'
    test_args.annotation_gtf = 'test/files/vep_test_files/annotation.gtf'
    test_args.output_prefix = 'test/files/vep_test_files/vep'
    test_args.verbose = True
    parse_vep(args=test_args)
