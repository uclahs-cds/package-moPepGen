""" A brute forth algorithm for calling variant peptides from a GVF file. """
import sys
import argparse
from typing import List, Dict, Tuple
from pathlib import Path
from itertools import combinations
from moPepGen import seqvar, aa, gtf, dna
from moPepGen.SeqFeature import FeatureLocation
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
        help='GVF file'
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
        miscleavage=args.miscleavage, min_mw=args.min_mw
    )

    variant_pool = seqvar.VariantRecordPool()
    with open(args.input_gvf) as handle:
        variant_pool.load_variants(handle, anno, genome)

    tx_id = list(variant_pool.transcriptional.keys())[0]

    tx_model = anno.transcripts[tx_id]
    tx_seq = tx_model.get_transcript_sequence(genome['chr1'])

    tx_peptides = tx_seq[tx_seq.orf.start:].translate(to_stop=True)\
        .enzymatic_cleave('trypsin', 'trypsin_exception')
    canonical_peptides = {str(peptide.seq) for peptide in tx_peptides}
    if tx_peptides and tx_peptides[0].seq.startswith('M'):
        canonical_peptides.add(str(tx_peptides[0].seq[1:]))
    variant_peptides = set()

    start_codon = FeatureLocation(
        start=tx_seq.orf.start, end=tx_seq.orf.start + 3
    )
    variants = variant_pool.transcriptional[tx_id]
    variants = [x for x in variants if x.location.start >= 3 and
        x.location < start_codon and not x.location.overlaps(start_codon)]

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

            has_start_altering = any(x.location.overlaps(start_codon) for x in comb)
            if not tx_model.is_cds_start_nf() and has_start_altering:
                cds_start = seq[tx_seq.orf.start + 1:].find('ATG')
                if cds_start == -1:
                    continue
            else:
                cds_start = tx_seq.orf.start
            aa_seq = seq[cds_start:].translate(to_stop=True)
            aa_seq = aa.AminoAcidSeqRecord(seq=aa_seq)
            peptides = aa_seq.enzymatic_cleave('trypsin', 'trypsin_exception')
            for peptide in peptides:
                if peptide is peptides[0] and peptide.seq.startswith('M'):
                    if str(peptide.seq[1:]) not in canonical_peptides:
                        variant_peptides.add(str(peptide.seq[1:]))
                if str(peptide.seq) not in canonical_peptides:
                    variant_peptides.add(str(peptide.seq))

    variant_peptides = list(variant_peptides)
    variant_peptides.sort()
    for peptide in variant_peptides:
        print(peptide, file=sys.stdout)
