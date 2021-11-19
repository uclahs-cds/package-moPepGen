""" A brute forth algorithm for calling variant peptides from a GVF file. """
import sys
import argparse
from typing import List, Dict, Tuple
from pathlib import Path
from itertools import combinations
from Bio import SeqUtils
from moPepGen import seqvar, aa, gtf, dna
from moPepGen.cli.common import add_args_cleavage, print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_brute_force(subparsers:argparse._SubParsersAction):
    """ parse command line arguments """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='bruteForce',
        help='Call variant peptide with the brute force algorithm.'
    )
    parser.add_argument(
        '-i', '--input-gvf',
        type=Path,
        help='GVF file',
        nargs="+"
    )
    parser.add_argument(
        '-r', '--reference-dir',
        type=Path,
        help='Reference directory. Must contain genome.fa, annotation.gtf, and'
        ' proteome.fasta. The directory should be generated by'
        ' the downsampleReference command'
    )
    add_args_cleavage(parser)
    parser.set_defaults(func=brute_force)
    print_help_if_missing_args(parser)
    return parser

def _parse_exclusion(exclusion) -> Tuple:
    """ Parse exclusion to values """
    start, ref, alt = exclusion.split('-')
    return int(start), ref, alt

def parse_variant_exclusion(exclusions:List[str]) -> Dict[Tuple,List[Tuple]]:
    """ Parse exclusion variants into """
    groups = {}
    for exclusion in exclusions:
        variant, targets = exclusion.split(':')
        targets = targets.split(',')
        variant = _parse_exclusion(variant)
        targets = [_parse_exclusion(target) for target in targets]
        groups[variant] = targets
    return groups

def brute_force(args):
    """ main """

    anno = gtf.GenomicAnnotation()
    anno.dump_gtf(args.reference_dir/'annotation.gtf')

    genome = dna.DNASeqDict()
    genome.dump_fasta(args.reference_dir/'genome.fasta')

    proteome = aa.AminoAcidSeqDict()
    proteome.dump_fasta(args.reference_dir/'proteome.fasta')

    rule = args.cleavage_rule
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    canonical_peptides = proteome.create_unique_peptide_pool(
        anno=anno, rule=rule, exception=exception,
        miscleavage=args.miscleavage, min_mw=args.min_mw,
        min_length=args.min_length, max_length=args.max_length
    )

    variant_pool = seqvar.VariantRecordPool()
    for gvf_file in args.input_gvf:
        with open(gvf_file) as handle:
            variant_pool.load_variants(handle, anno, genome)

    tx_id = list(variant_pool.transcriptional.keys())[0]

    tx_model = anno.transcripts[tx_id]
    tx_seq = tx_model.get_transcript_sequence(genome['chr1'])

    variant_peptides = set()

    start_index = tx_seq.orf.start + 3
    variants = []

    for variant in variant_pool.transcriptional[tx_id]:
        if variant.location.start < start_index -1:
            continue
        if variant.location.start == start_index - 1:
            variant.to_end_inclusion(tx_seq)
        variants.append(variant)

    for i in range(len(variants)):
        for comb in combinations(variants, i + 1):
            skip = False

            for j, left in enumerate(comb):
                if j == len(comb) - 1:
                    continue
                for right in comb[j+1:]:
                    if left.location.end >= right.location.start:
                        skip = True

            if skip:
                continue

            seq = tx_seq.seq
            offset = 0
            for variant in comb:
                start = variant.location.start + offset
                end = variant.location.end + offset
                offset = offset + len(variant.alt) - len(variant.ref)
                seq = seq[:start] + variant.alt + seq[end:]

            if not tx_seq.orf:
                cds_start_positions = tx_seq.find_all_start_codons()
            else:
                cds_start = tx_seq.orf.start
                cds_start_positions = [cds_start]
            for cds_start in cds_start_positions:
                cds_stop = len(seq) - (len(seq) - cds_start) % 3
                aa_seq = seq[cds_start:cds_stop].translate(to_stop=True)
                aa_seq = aa.AminoAcidSeqRecord(seq=aa_seq)
                peptides = aa_seq.enzymatic_cleave('trypsin', 'trypsin_exception')
                for peptide in peptides:
                    if peptide is peptides[0] and peptide.seq.startswith('M'):
                        is_valid = peptide_is_valid(peptide.seq[1:],
                            args.min_length, args.max_length, args.min_mw)
                        if is_valid and \
                                str(peptide.seq[1:]) not in canonical_peptides:
                            variant_peptides.add(str(peptide.seq[1:]))
                    is_valid = peptide_is_valid(peptide.seq, args.min_length,
                        args.max_length, args.min_mw)
                    if is_valid and str(peptide.seq) not in canonical_peptides:
                        variant_peptides.add(str(peptide.seq))

    variant_peptides = list(variant_peptides)
    variant_peptides.sort()
    for peptide in variant_peptides:
        print(peptide, file=sys.stdout)

def peptide_is_valid(peptide:str, min_length:int, max_length:int, min_mw:float
        ) -> bool:
    """ Check whether the peptide is valid """
    return min_length <= len(peptide) <= max_length and \
        SeqUtils.molecular_weight(peptide, 'protein') >= min_mw
