""" A brute forth algorithm for calling variant peptides from a GVF file. """
import sys
import argparse
from typing import Iterable, List, Dict, Tuple
from pathlib import Path
from itertools import combinations
from Bio import SeqUtils
from moPepGen import gtf, seqvar, aa, dna, params
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.cli.common import add_args_cleavage, print_help_if_missing_args
from moPepGen.util.common import load_references


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
    parser.add_argument(
        '-f', '--force',
        action='store_true',
        help='If not set, the program stops when there are more than 10'
        ' variants. When this flag is set, the program runs anyway. Noted that '
        ' the runtime is going to increase quickly after 10 variants.',
        default=False
    )
    parser.add_argument(
        '--variant-ids',
        type=str,
        help='List of variant labels.',
        nargs='*'
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


class BruteForceVariantPeptideCaller():
    """ Variant peptide caller using the brute force algorithm. """
    def __init__(self, reference_data:params.ReferenceData=None,
            cleavage_params:params.CleavageParams=None,
            variants:List[seqvar.VariantRecord]=None,
            canonical_peptides=None, tx_id:str=None,
            tx_model:gtf.TranscriptAnnotationModel=None,
            tx_seq:dna.DNASeqRecord=None, start_index:int=None):
        """ Constructor """
        self.reference_data = reference_data
        self.cleavage_params = cleavage_params
        self.variants = variants or []
        self.canonical_peptides = canonical_peptides
        self.tx_id = tx_id
        self.tx_model = tx_model
        self.tx_seq = tx_seq
        self.start_index = start_index

    def create_canonical_peptide_pool(self):
        """ Create canonical peptide pool. """
        proteome = self.reference_data.proteome
        par = self.cleavage_params
        self.canonical_peptides = proteome.create_unique_peptide_pool(
            anno=self.reference_data.anno,
            rule=par.enzyme,
            exception=par.exception,
            miscleavage=par.miscleavage,
            min_mw=par.min_mw,
            min_length=par.min_length,
            max_length=par.max_length
        )

    def load_variant_records(self, variant_pool:seqvar.VariantRecordPool,
            variant_ids:Iterable[str]=None):
        """ Load variant records associated with the particular transcript. """
        series = variant_pool[self.tx_id]
        variants = []
        for variant in series.transcriptional:
            if variant.location.start < self.start_index -1:
                continue
            if self.tx_model.is_mrna_end_nf() and variant.location.end <= self.tx_seq.orf.end - 3:
                continue
            if variant_ids:
                if variant.id not in variant_ids:
                    continue
            if variant.location.start == self.start_index - 1:
                variant.to_end_inclusion(self.tx_seq)
            variants.append(variant)
        self.variants = variants

    def get_start_index(self):
        """ Get the "start index" used for filtering variants. For noncoding
        transcript,s the  "start index" is set to the very beginning. """
        if self.tx_seq.orf:
            self.start_index = self.tx_seq.orf.start + 3
        else:
            self.start_index = 3

    def peptide_is_valid(self, peptide:str) -> bool:
        """ Check whether the peptide is valid """
        if self.canonical_peptides and peptide in self.canonical_peptides:
            return False
        min_len = self.cleavage_params.min_length
        max_len = self.cleavage_params.max_length
        min_mw = self.cleavage_params.min_mw
        return min_len <= len(peptide) <= max_len \
            and SeqUtils.molecular_weight(peptide, 'protein') >= min_mw

    def has_any_variant(self, lhs:int, rhs:int, cds_start:int,
            variants:List[seqvar.VariantRecord]) -> bool:
        """ Check whether the given range of the transcript has any variant
        associated. """
        offset = 0
        query = FeatureLocation(start=lhs, end=rhs)
        start_loc = FeatureLocation(start=cds_start, end=cds_start + 3)
        for variant in variants:
            var_size = len(variant.alt) - len(variant.ref)
            loc = FeatureLocation(
                start=variant.location.start + offset,
                end=variant.location.end + offset + var_size
            )
            if loc.start > rhs + 3:
                break
            is_start_gain = start_loc.overlaps(loc)
            is_frameshifting = cds_start < loc.start and variant.is_frameshifting()
            is_cleavage_gain = loc.overlaps(FeatureLocation(start=lhs - 3, end=lhs)) \
                    or loc.overlaps(FeatureLocation(start=rhs, end=rhs + 3)) \
                if cds_start != lhs \
                else loc.overlaps(FeatureLocation(start=rhs, end=rhs + 3))
            orf_index = cds_start % 3
            codon_pos = loc.start - (loc.start - orf_index) % 3
            is_stop_lost = cds_start < loc.start \
                    and self.tx_seq.seq[codon_pos:codon_pos + 3] in ['TAA', 'TAG', 'TGA']
            if loc.overlaps(query) \
                    or is_start_gain \
                    or is_frameshifting \
                    or is_cleavage_gain \
                    or is_stop_lost:
                return True
            offset += var_size
        return False

    def have_incompatible_variants(self, variants:Iterable[seqvar.VariantRecord]) -> bool:
        """ Whether there are incompatible variants, i.e. variants that overlap
        with each other. """
        for j, left in enumerate(variants):
            if j == len(variants) - 1:
                continue
            for right in variants[j+1:]:
                if left.location.end >= right.location.start:
                    return True
        if self.tx_model.is_mrna_end_nf():
            for variant in variants:
                if variant.location.end >= self.tx_seq.orf.end:
                    return True
        return False

    def call_peptides(self, force:bool):
        """ Call variant peptides """
        if not self.variants:
            return

        variant_peptides = set()

        if len(self.variants) > 10 and not force:
            raise ValueError(
                f"{len(self.variants)} variants is too many for this brute force algorithm."
            )

        variants = self.variants
        tx_model = self.tx_model
        tx_seq = self.tx_seq
        rule = self.cleavage_params.enzyme
        exception = self.cleavage_params.exception
        for i in range(len(self.variants)):
            for comb in combinations(variants, i + 1):
                if self.have_incompatible_variants(comb):
                    continue

                seq = tx_seq.seq
                offset = 0
                if tx_model.is_protein_coding and tx_model.is_mrna_end_nf():
                    cur_cds_end = tx_seq.orf.end

                for variant in comb:
                    start = variant.location.start + offset
                    end = variant.location.end + offset
                    if tx_model.is_protein_coding and tx_model.is_mrna_end_nf():
                        cur_cds_end += len(variant.alt) - len(variant.ref)
                    offset = offset + len(variant.alt) - len(variant.ref)
                    seq = seq[:start] + variant.alt + seq[end:]

                if not (tx_model.is_protein_coding and tx_model.is_mrna_end_nf()):
                    cur_cds_end = len(seq)

                if not tx_model.is_protein_coding:
                    alt_seq = dna.DNASeqRecord(seq)
                    cds_start_positions = alt_seq.find_all_start_codons()
                else:
                    cds_start = tx_seq.orf.start
                    cds_start_positions = [cds_start]

                for cds_start in cds_start_positions:
                    if not tx_model.is_protein_coding:
                        cur_cds_end = len(seq)
                    else:
                        cur_cds_end = cur_cds_end - (cur_cds_end - cds_start) % 3

                    aa_seq = seq[cds_start:cur_cds_end].translate(to_stop=True)
                    aa_seq = aa.AminoAcidSeqRecord(seq=aa_seq)

                    sites = aa_seq.find_all_enzymatic_cleave_sites(rule, exception)
                    sites.insert(0, 0)
                    sites.append(len(aa_seq))
                    for j, lhs in enumerate(sites[:-1]):
                        for k in range(j + 1, min([j + 3, len(sites) - 1]) + 1):
                            rhs = sites[k]
                            tx_lhs = cds_start + lhs * 3
                            tx_rhs = cds_start + rhs * 3
                            if not self.has_any_variant(tx_lhs, tx_rhs, cds_start, comb):
                                continue

                            peptides = [aa_seq[lhs:rhs]]
                            if peptides[0].seq.startswith('M') and lhs == 0:
                                peptides.append(peptides[0][1:])
                            for peptide in peptides:
                                is_valid = self.peptide_is_valid(peptide.seq)
                                if is_valid:
                                    variant_peptides.add(str(peptide.seq))

        variant_peptides = list(variant_peptides)
        variant_peptides.sort()
        for peptide in variant_peptides:
            print(peptide, file=sys.stdout)

def brute_force(args):
    """ main """
    # Load genomic references
    anno, genome, proteome = load_references(
        path_anno=args.reference_dir/'annotation.gtf',
        path_genome=args.reference_dir/'genome.fasta',
        path_proteome=args.reference_dir/'proteome.fasta'
    )
    reference_data = params.ReferenceData(
        genome=genome,
        anno=anno,
        proteome=proteome,
        canonical_peptides=None
    )

    # load GVF files
    variant_pool = seqvar.VariantRecordPool()
    for path in args.input_gvf:
        with open(path) as handle:
            variant_pool.load_variants(
                handle=handle,
                anno=reference_data.anno,
                genome=reference_data.genome
            )
    tx_id = list(variant_pool.data.keys())[0]

    caller = BruteForceVariantPeptideCaller()
    caller.reference_data = reference_data
    caller.tx_id = tx_id
    caller.tx_model = caller.reference_data.anno.transcripts[caller.tx_id]
    caller.tx_seq = caller.tx_model.get_transcript_sequence(caller.reference_data.genome['chr1'])

    caller.cleavage_params = params.CleavageParams(
        enzyme=args.cleavage_rule,
        exception='trypsin_exception' if args.cleavage_rule == 'trypsin' else None,
        miscleavage=int(args.miscleavage),
        min_mw=float(args.min_mw),
        min_length=args.min_length,
        max_length=args.max_length
    )

    caller.create_canonical_peptide_pool()
    caller.get_start_index()
    caller.load_variant_records(variant_pool)

    caller.call_peptides(force=args.force)
