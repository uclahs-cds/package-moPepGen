""" A brute forth algorithm for calling variant peptides from a TVF file. """
import sys
import argparse
from typing import List, Dict, Tuple
from pathlib import Path
from itertools import combinations
from moPepGen import seqvar, aa, gtf, dna


def parse_args():
    """ parse command line argments """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input-tvf',
        type=Path,
        help='TVF file'
    )
    parser.add_argument(
        '-r', '--reference-dir',
        type=Path,
        help='Index directory'
    )
    parser.add_argument(
        '-c', '--cleavage-rule',
        type=str,
        help='Cleavage rule. Defaults to trypsin.',
        default='trypsin',
        metavar=''
    )
    parser.add_argument(
        '-m', '--miscleavage',
        type=int,
        help='Number of cleavages to allow. Defaults to 2.',
        metavar='',
        default=2
    )
    parser.add_argument(
        '-w', '--min-mw',
        type=float,
        help='The minimal molecular weight of the non-canonical peptides.'
        'Defaults to 500',
        default=500.,
        metavar=''
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

def main(args):
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
        miscleavage=args.miscleavage, min_mw=args.min_mw
    )

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

            for j, left in enumerate(comb):
                if j == len(comb) - 1:
                    continue
                for right in comb[j+1:]:
                    if left.location.end > right.location.start:
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
                if peptide is peptides[0] and not tx_model.is_cds_start_nf() \
                        and peptide.seq.startswith('M'):
                    if str(peptide.seq[1:]) not in canonical_peptides:
                        variant_peptides.add(str(peptide.seq[1:]))
                if str(peptide.seq) not in canonical_peptides:
                    variant_peptides.add(str(peptide.seq))

    variant_peptides = list(variant_peptides)
    variant_peptides.sort()
    for peptide in variant_peptides:
        print(peptide, file=sys.stdout)


if __name__ == '__main__':
    _args = parse_args()
    main(_args)
