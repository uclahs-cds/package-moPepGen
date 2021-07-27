import sys
import argparse
from typing import List, Dict, Tuple
from pathlib import Path
import pickle
from moPepGen import seqvar, aa
from itertools import combinations


def parse_args():
    """ parse command line argments """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input-tvf',
        type=Path,
        help='TVF file'
    )
    parser.add_argument(
        '-d', '--index-dir',
        type=Path,
        help='Index directory'
    )
    parser.add_argument(
        '-x', '--exclusion',
        type=str,
        help='Variant exclusions. Format:'
        '<variant1>:<exclusion-1>,<exclusion-2> '
        '<variant2>:<exclusion-3>,<exclusion-4>. Variants are represented as '
        '<start>-<ref>-<alt>',
        nargs="*"
    )
    return parser.parse_args()

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

def tuplify_variant(variant:seqvar.VariantRecord) -> str:
    """ Convert a variant to a tuple """
    return int(variant.location.start), variant.ref, variant.alt

def main(args):
    """ main """
    if args.exclusion:
        exclusion = parse_variant_exclusion(args.exclusion)
    else:
        exclusion = {}

    with open(args.index_dir/'genome.pickle', 'rb') as handle:
        genome = pickle.load(handle)
    with open(args.index_dir/'annotation.pickle', 'rb') as handle:
        anno = pickle.load(handle)

    variants = list(seqvar.io.parse(args.input_tvf))
    variants.sort()
    tx_id = variants[0].location.seqname
    tx_model = anno.transcripts[tx_id]
    tx_seq = tx_model.get_transcript_sequence(genome['chr1'])

    canonical_peptides = tx_seq[tx_seq.orf.start:].translate(to_stop=True)\
        .enzymatic_cleave('trypsin', 'trypsin_exception')
    canonical_peptides = {str(peptide.seq) for peptide in canonical_peptides}
    variant_peptides = set()

    for i in range(len(variants)):
        for comb in combinations(variants, i + 1):
            skip = False
            for variant in comb:
                target = tuplify_variant(variant)
                if target not in exclusion:
                    continue
                exclude_variants = exclusion[target]
                for x in variants:
                    if tuplify_variant(x) in exclude_variants:
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
            aa_seq = seq[tx_seq.orf.start:].translate(to_stop=True)
            aa_seq = aa.AminoAcidSeqRecord(seq=aa_seq)
            peptides = aa_seq.enzymatic_cleave('trypsin', 'trypsin_exception')
            for peptide in peptides:
                if str(peptide.seq) not in canonical_peptides:
                    variant_peptides.add(str(peptide.seq))

    variant_peptides = list(variant_peptides)
    variant_peptides.sort()
    for peptide in variant_peptides:
        print(peptide, file=sys.stdout)


if __name__ == '__main__':
    args = argparse.Namespace()
    args.input_tvf = Path('test/files/vep/CPCG0103_gencode_aa_indel_ENST00000314675.11.tvf')
    args.index_dir = Path('test/files/downsampled_index/ENST00000314675.11')
    args.exclusion = ['1213-TCCCGGCCCCA-T:1204-CGGTCCCGG-C,1219-CCCCAGTCCCGGTCCCGG-C']
    # args = parse_args()
    main(args)

