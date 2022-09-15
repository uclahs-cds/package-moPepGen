""" A brute forth algorithm for calling variant peptides from a GVF file. """
import sys
import argparse
from typing import Iterable, List, Dict, Set, Tuple
from pathlib import Path
from itertools import combinations
from Bio import SeqUtils
from Bio.Seq import Seq
from moPepGen import gtf, seqvar, aa, dna, params
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.cli.common import add_args_cleavage, print_help_if_missing_args
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool
from moPepGen.seqvar.VariantRecordWithCoordinate import VariantRecordWithCoordinate
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
            variant_pool:seqvar.VariantRecordPool=None,
            canonical_peptides=None, tx_id:str=None,
            tx_model:gtf.TranscriptAnnotationModel=None,
            tx_seq:dna.DNASeqRecord=None, gene_seq:dna.DNASeqRecord=None,
            start_index:int=None, variant_peptides:Set[str]=None):
        """ Constructor """
        self.reference_data = reference_data
        self.cleavage_params = cleavage_params
        self.variant_pool = variant_pool or VariantRecordPool()
        self.canonical_peptides = canonical_peptides
        self.tx_id = tx_id
        self.tx_model = tx_model
        self.tx_seq = tx_seq
        self.gene_seq = gene_seq
        self.start_index = start_index
        self.variant_peptides = variant_peptides or set()

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

    def get_gene_seq(self) -> dna.DNASeqRecord:
        """ Get the gene sequence and cache it if it is not already cached. """
        if self.gene_seq:
            return self.gene_seq
        gene_id = self.tx_model.gene_id
        gene_model = self.reference_data.anno.genes[gene_id]
        chrom = gene_model.chrom
        self.gene_seq = gene_model.get_gene_sequence(self.reference_data.genome[chrom])
        return self.gene_seq

    def load_relevant_variants(self, pool:VariantRecordPool):
        """ Load relevant variants. Relevant variants are those associated with
        the transcript of `self.tx_id`, or the donor transcript of any fusion. """
        if self.tx_id not in pool:
            return

        for variant in pool[self.tx_id].transcriptional:
            self.variant_pool.add_transcriptional_variant(variant)
        for variant in pool[self.tx_id].intronic:
            self.variant_pool.add_intronic_variant(variant)
        for variant in pool[self.tx_id].fusion:
            self.variant_pool.add_fusion_variant(variant)
            right_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
            if right_tx_id not in pool:
                continue
            for v in pool[right_tx_id].transcriptional:
                self.variant_pool.add_transcriptional_variant(v)
            for v in pool[right_tx_id].intronic:
                self.variant_pool.add_intronic_variant(v)
        for variant in pool[self.tx_id].circ_rna:
            self.variant_pool.add_circ_rna(variant)

    def get_variants(self, tx_id:str, start:int, end:int,
            variant_ids:Iterable[str]=None) -> List[seqvar.VariantRecord]:
        """ Load variant records associated with the particular transcript. """
        series = self.variant_pool[tx_id]
        variants = []
        for variant in series.transcriptional:
            if variant.location.start < start -1:
                continue
            if self.tx_model.is_mrna_end_nf() and variant.location.end <= end:
                continue
            if variant_ids:
                if variant.id not in variant_ids:
                    continue
            if variant.location.start == start - 1:
                variant.to_end_inclusion(self.tx_seq)
            variants.append(variant)
        return variants

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

    def is_stop_lost(self, variant:seqvar.VariantRecord, reading_frame_index:int) -> bool:
        """ Check whether the variant is a stop lost mutation. """
        orf_index = self.start_index % 3

        if self.tx_model.is_protein_coding:
            if orf_index != reading_frame_index:
                return False
            orf_end = self.tx_seq.orf.end
            stop_codon = FeatureLocation(orf_end, orf_end + 3)
            return variant.location.overlaps(stop_codon)

        i = variant.variant.location.start - (variant.location.start - orf_index) % 3
        while i < variant.location.end:
            if self.tx_seq.seq[i:i+3] in ['TAA', 'TAG', 'TGA']:
                return True
            i += 3
        return False

    @staticmethod
    def has_any_variant(lhs:int, rhs:int, cds_start:int, prev_cds_start:int,
            variants:List[seqvar.VariantRecordWithCoordinate],
            variants_stop_lost:List[Tuple[bool,bool,bool]],
            variants_silent_mutation:List[Tuple[bool,bool,bool]]) -> bool:
        """ Check whether the given range of the transcript has any variant
        associated. """
        offset = 0
        query = FeatureLocation(start=lhs, end=rhs)
        start_loc = FeatureLocation(start=cds_start, end=cds_start + 3)
        for i, variant_coordinate in enumerate(variants):
            variant = variant_coordinate.variant
            loc = variant_coordinate.location
            if loc.start > rhs + 3:
                break
            is_start_gain = start_loc.overlaps(loc)
            is_frameshifting = prev_cds_start < loc.start and variant.is_frameshifting()
            is_cleavage_gain = loc.overlaps(FeatureLocation(start=lhs - 3, end=lhs)) \
                    or loc.overlaps(FeatureLocation(start=rhs, end=rhs + 3)) \
                if cds_start != lhs \
                else loc.overlaps(FeatureLocation(start=rhs, end=rhs + 3))

            is_stop_lost = variants_stop_lost[i][cds_start % 3] \
                and variants[i].location.start > cds_start

            is_silent_mutation = variants_silent_mutation[i][cds_start % 3] \
                and variants[i].location.start > cds_start

            if (loc.overlaps(query) \
                        or is_start_gain \
                        or is_frameshifting \
                        or is_cleavage_gain \
                        or is_stop_lost )\
                    and not is_silent_mutation:
                return True
            offset += len(variant.alt) - len(variant.ref)
        return False

    @staticmethod
    def has_overlapping_variants(variants:List[VariantRecord]) -> bool:
        """ Checks if any variants overlap. """
        for i, left in enumerate(variants):
            if i == len(variants) - 1:
                continue
            for right in variants[i+1:]:
                if left.location.end >= right.location.start:
                    return True
        return False

    def has_any_invalid_variants_on_inserted_sequences(self,
            pool:seqvar.VariantRecordPool) -> bool:
        """ Checks if any variants carried by fusion or alt splicing. In valid
        variants are those not in the region of sequence introduced by fusion.
        For example, a variant of the donor transcript before the breakpoint is
        invalid.
        """
        alt_splices = [x for x in pool[self.tx_id].transcriptional if x.is_alternative_splicing()]
        if not pool[self.tx_id].fusion and not alt_splices:
            return False

        inserted_intronic_region:Dict[str, List[FeatureLocation]] = {}

        if pool[self.tx_id].fusion:
            fusion = pool[self.tx_id].fusion[0]
            left_insert_start = fusion.attrs['LEFT_INSERTION_START']
            left_insert_end = fusion.attrs['LEFT_INSERTION_END']
            right_insert_start = fusion.attrs['RIGHT_INSERTION_START']
            right_insert_end = fusion.attrs['RIGHT_INSERTION_END']
            right_tx_id = fusion.attrs['ACCEPTER_TRANSCRIPT_ID']
            right_breakpoint = fusion.get_accepter_position()

            if right_tx_id in pool:
                alt_splices += [x for x in pool[right_tx_id].transcriptional
                    if x.is_alternative_splicing()]

            if any(x.location.start > fusion.location.start
                    for x in pool[self.tx_id].transcriptional):
                return True

            if left_insert_start:
                loc = FeatureLocation(start=left_insert_start, end=left_insert_end)
                if self.tx_id not in inserted_intronic_region:
                    inserted_intronic_region[self.tx_id] = []
                inserted_intronic_region[self.tx_id].append(loc)

            if right_insert_start:
                loc = FeatureLocation(start=right_insert_start, end=right_insert_end)
                if self.tx_id not in inserted_intronic_region:
                    inserted_intronic_region[self.tx_id] = []
                inserted_intronic_region[self.tx_id].append(loc)

            if right_tx_id in pool \
                    and any(x.location.start < right_breakpoint
                        for x in pool[right_tx_id].transcriptional):
                return True

        for alt_splice in alt_splices:
            donor_start = alt_splice.attrs.get('DONOR_START')
            if not donor_start:
                continue
            donor_start = int(donor_start)
            donor_end = int(alt_splice.attrs['DONOR_END'])
            loc = FeatureLocation(start=donor_start, end=donor_end)
            if alt_splice.transcript_id not in inserted_intronic_region:
                inserted_intronic_region[alt_splice.transcript_id] = []
            inserted_intronic_region[alt_splice.transcript_id].append(loc)

        for tx_id in pool:
            if tx_id not in inserted_intronic_region:
                return True
            for v in pool[tx_id].intronic:
                if not any(x.is_superset(v.location) for x in inserted_intronic_region[tx_id]):
                    return True
        return False

    @staticmethod
    def has_any_invalid_variant_on_circ() -> bool:
        """ Checks if any variants are invlid with circRNA. """
        return False

    def has_incompatible_variants(self, pool:seqvar.VariantRecordPool) -> bool:
        """ Whether there are incompatible variants, i.e. variants that overlap
        with each other. """
        # check if there is any variants
        if self.tx_id not in pool:
            return True

        if not pool[self.tx_id].transcriptional\
                and not pool[self.tx_id].circ_rna \
                and not pool[self.tx_id].fusion:
            return True

        # check if there are multiple novel transcript variants (fusion + circ)
        if len(pool[self.tx_id].circ_rna) + len(pool[self.tx_id].circ_rna) > 1:
            return True

        if self.has_any_invalid_variants_on_inserted_sequences(pool):
            return True

        # if self.has_any_invalid_variant_on_circ(pool):
        #     return True

        for series in pool.data.values():
            if self.has_overlapping_variants(series.transcriptional):
                return True
            if self.has_overlapping_variants(series.intronic):
                return True

        if self.tx_model.is_mrna_end_nf():
            for variant in pool[self.tx_id].transcriptional:
                if variant.location.end >= self.tx_seq.orf.end:
                    return True
        return False

    @staticmethod
    def find_prev_cds_start_same_frame(cds_start:int, cds_start_positions:List[int]):
        """ find he previous cds start site in the same reading frame. """
        if cds_start == 0:
            return -1
        reading_frame_index = cds_start % 3
        for site in cds_start_positions[::-1]:
            if site >= cds_start:
                continue
            if site % 3 == reading_frame_index:
                return site
        return -1

    def get_variant_sequence(self, seq:Seq, location:FeatureLocation,
            offset:int, variants:List[seqvar.VariantRecord],
            pool:seqvar.VariantRecordPool
            ) -> Tuple[Seq, List[seqvar.VariantRecordWithCoordinate]]:
        """ Get variant sequence. """
        var_seq = seq
        variant_coordinates = []
        local_offset = 0
        for variant in variants:
            if variant.is_alternative_splicing():
                if variant.type == 'Deletion':
                    start = variant.location.start + local_offset - location.start
                    end = variant.location.end + local_offset - location.start
                    loc = FeatureLocation(start=start, end=start + 1)
                    variant_coordinate = seqvar.VariantRecordWithCoordinate(
                        variant=variant,
                        location=loc
                    )

                    variant_coordinates.append(variant_coordinate)
                    alt_seq = var_seq[start:start+1]
                    local_offset = local_offset + len(alt_seq) - len(variant.location)
                    var_seq = var_seq[:start] + alt_seq + var_seq[end:]

                elif variant.type == 'Insertion':
                    start = variant.location.start + local_offset - location.start
                    end = variant.location.end + local_offset - location.start

                    gene_seq = self.get_gene_seq()
                    donor_start = variant.get_donor_start()
                    donor_end = variant.get_donor_end()
                    alt_seq = str(gene_seq.seq[donor_start:donor_end])
                    loc = FeatureLocation(start=donor_start, end=donor_end)
                    insert_variants = [x for x in pool[self.tx_id].intronic
                        if loc.is_superset(x.location)]
                    alt_seq, insert_variants = self.get_variant_sequence(
                        seq=alt_seq, location=loc, offset=start,
                        variants=insert_variants, pool=pool
                    )

                    variant_coordinate = seqvar.VariantRecordWithCoordinate(
                        variant=variant,
                        location=FeatureLocation(start=start, end=start + len(alt_seq))
                    )
                    variant_coordinates.append(variant_coordinate)
                    variant_coordinates += insert_variants
                    local_offset = local_offset + len(alt_seq) + 1 - len(variant.location)
                    var_seq = var_seq[:start+1] + alt_seq + var_seq[end:]

                elif variant.type == 'Substitution':
                    start = variant.location.start + local_offset - location.start
                    end = variant.location.end + local_offset - location.start

                    gene_seq = self.get_gene_seq()
                    donor_start = variant.get_donor_start()
                    donor_end = variant.get_donor_end()
                    alt_seq = str(gene_seq.seq[donor_start:donor_end])
                    loc = FeatureLocation(start=donor_start, end=donor_end)
                    insert_variants = [x for x in pool[self.tx_id].intronic
                        if loc.is_superset(x.location)]
                    alt_seq, insert_variants = self.get_variant_sequence(
                        seq=alt_seq, location=loc, offset=start,
                        variants=insert_variants, pool=pool
                    )

                    variant_coordinate = seqvar.VariantRecordWithCoordinate(
                        variant=variant,
                        location=FeatureLocation(start=start, end=start + len(alt_seq))
                    )
                    variant_coordinates.append(variant_coordinate)
                    variant_coordinates += insert_variants
                    local_offset = local_offset + len(alt_seq) - len(variant.location)
                    var_seq = var_seq[:start] + alt_seq + var_seq[end:]

            else:
                start = variant.location.start + local_offset - location.start
                end = variant.location.end + local_offset - location.start
                loc = FeatureLocation(
                    start=start + offset, end=start + len(variant.alt) + offset
                )
                variant_coordinate = seqvar.VariantRecordWithCoordinate(
                    variant=variant,
                    location=loc
                )

                variant_coordinates.append(variant_coordinate)
                local_offset = local_offset + len(variant.alt) - len(variant.ref)
                var_seq = var_seq[:start] + variant.alt + var_seq[end:]

        return var_seq, variant_coordinates

    def get_variant_sequence_fusion(self, seq:Seq, variants:seqvar.VariantRecordPool
            ) -> Tuple[Seq, List[seqvar.VariantRecordWithCoordinate]]:
        """ Get variant sequence with fusion. """
        number_of_fusion = len(variants[self.tx_id].fusion)
        if not number_of_fusion == 1:
            raise ValueError(
                f"Should have exactly 1 fusion, but {number_of_fusion} were found."
            )
        fusion = variants[self.tx_id].fusion[0]
        var_seq = seq[:fusion.location.start]
        location = FeatureLocation(start=0, end=len(var_seq))
        var_seq, variant_coordinates = self.get_variant_sequence(
            seq=var_seq, location=location, local_offset=len(var_seq),
            variants=variants[self.tx_id].transcriptional, pool=variants
        )

        left_insert_start = fusion.attrs['LEFT_INSERTION_START']
        left_insert_end = fusion.attrs['LEFT_INSERTION_END']
        right_insert_start = fusion.attrs['RIGHT_INSERTION_START']
        right_insert_end = fusion.attrs['RIGHT_INSERTION_END']
        right_tx_id = fusion.attrs['ACCEPTER_TRANSCRIPT_ID']
        right_gene_id = fusion.attrs['ACCEPTER_GENE_ID']

        additional_seq = Seq('')
        additional_variants:List[VariantRecordWithCoordinate] = []

        # left insertion
        if left_insert_start is not None:
            gene_seq = self.get_gene_seq()
            location = FeatureLocation(start=left_insert_start, end=left_insert_end)
            insert_seq = gene_seq.seq[left_insert_start:left_insert_end]
            insert_variants = [x for x in variants[self.tx_id].intronic
                if location.is_superset(x.location)]
            insert_seq, insert_variants = self.get_variant_sequence(
                seq=insert_seq, location=location,
                local_offset=len(var_seq) + len(additional_seq),
                variants=insert_variants, pool=variants
            )
            additional_seq += insert_seq
            additional_variants += insert_variants

        # right insertion
        if right_insert_start is not None:
            gene_model = self.reference_data.anno.genes[right_gene_id]
            chrom = gene_model.chrom
            gene_seq = gene_model.get_gene_sequence(self.reference_data.genome[chrom])
            location = FeatureLocation(start=right_insert_start, end=right_insert_end)
            insert_seq = gene_seq.seq[right_insert_start:right_insert_end]
            insert_seq, insert_variants = self.get_variant_sequence(
                seq=insert_seq, location=location,
                local_offset=len(var_seq) + len(additional_seq),
                variants=variants[right_tx_id].intronic if right_tx_id in variants else [],
                pool=variants
            )
            additional_seq += insert_seq
            additional_variants += insert_variants

        right_tx_model = self.reference_data.anno.transcripts[right_tx_id]
        accepter_chrom = right_tx_model.transcript.location.seqname
        breakpoint_gene = fusion.get_accepter_position()
        breakpoint_tx = self.reference_data.anno.coordinate_gene_to_transcript(
            index=breakpoint_gene,
            gene=right_gene_id,
            transcript=right_tx_id
        )
        accepter_seq = right_tx_model.get_transcript_sequence(
            self.reference_data.genome[accepter_chrom]
        )
        location = FeatureLocation(start=breakpoint_tx, end=len(accepter_seq))
        insert_seq = accepter_seq.seq[breakpoint_tx:]
        insert_seq, insert_variants = self.get_variant_sequence(
            seq=insert_seq, location=location,
            local_offset=len(var_seq) + len(additional_seq),
            variants=variants[right_tx_id].transcriptional if right_tx_id in variants else [],
            pool=variants
        )
        additional_seq += insert_seq
        additional_variants += insert_variants

        location = FeatureLocation(start=len(var_seq), end=len(var_seq) + len(additional_seq))
        fusion_var = VariantRecordWithCoordinate(variant=fusion, location=location)
        var_seq += additional_seq
        variant_coordinates.append(fusion_var)
        variant_coordinates += additional_variants

        return var_seq, variant_coordinates

    def get_variant_ref_seq(self, variant:seqvar.VariantRecord) -> Seq:
        """ Get the reference sequence of a variant """
        if variant.type in ['Deletion', 'Substitution']:
            return self.tx_seq.seq[variant.location.start:variant.location.end]
        return variant.ref

    def check_stop_lost(self, seq:str, variants:List[seqvar.VariantRecordWithCoordinate]
            ) -> List[Tuple[bool, bool, bool]]:
        """ Checks if each variant is stop lost on each reading frame. """
        stop_lost:List[Tuple[bool,bool,bool]] = []
        for variant in variants:
            if variant.variant.is_fusion() \
                    or variant.variant.is_circ_rna() \
                    or variant.variant.transcript_id != self.tx_id \
                    or 'TRANSCRIPT_ID' in variant.variant.attrs:
                stop_lost.append((False, False, False))
                continue
            stop_lost_i:List[bool] = []
            for i in range(3):
                # lhs & rhs: position of the first (lhs) and last (lhs) codon
                # that overlaps with the variant
                lhs = variant.location.start - (variant.location.start - i) % 3
                var_ref = self.get_variant_ref_seq(variant.variant)
                ref_seq = seq[lhs:variant.location.start] + var_ref
                n_carry_over = 3 - (len(ref_seq) % 3)
                rhs = min(len(seq), variant.location.end + n_carry_over)
                ref_seq += seq[variant.location.end:rhs]
                j = 0
                while j + 3 <= len(ref_seq):
                    if ref_seq[j:j+3] in ['TAA', 'TAG', 'TGA']:
                        stop_lost_i.append(True)
                        break
                    j += 3
                else:
                    stop_lost_i.append(False)
            stop_lost.append(tuple(stop_lost_i))
        return stop_lost

    def check_silent_mutation(self, seq:str,
            variants:List[seqvar.VariantRecordWithCoordinate]):
        """ """
        silent_mutation:List[Tuple[bool,bool,bool]] = []
        for variant in variants:
            if variant.variant.is_fusion() \
                    or variant.variant.is_circ_rna() \
                    or variant.variant.is_alternative_splicing():
                silent_mutation.append((False, False, False))
                continue
            silent_i:List[bool] = []
            for i in range(3):
                # lhs & rhs: position of the first (lhs) and last (lhs) codon
                # that overlaps with the variant
                lhs = variant.location.start - (variant.location.start - i) % 3
                var_ref = self.get_variant_ref_seq(variant.variant)
                ref_seq = seq[lhs:variant.location.start] + var_ref
                n_carry_over = 3 - (len(ref_seq) % 3)
                rhs = min(len(seq), variant.location.end + n_carry_over)
                ref_seq += seq[variant.location.end:rhs]

                var_seq = seq[lhs:variant.location.end]
                n_carry_over = 3 - (len(var_seq) % 3)
                rhs = min(len(seq), variant.location.end + n_carry_over)
                var_seq += seq[variant.location.end:rhs]

                silent_i.append(
                    ref_seq.translate(to_stop=False) == var_seq.translate(to_stop=False)
                )
            silent_mutation.append(tuple(silent_i))
        return silent_mutation

    def call_peptides_main(self, variants:seqvar.VariantRecordPool):
        """ Call peptide main """
        tx_model = self.tx_model
        tx_seq = self.tx_seq
        rule = self.cleavage_params.enzyme
        exception = self.cleavage_params.exception

        if variants[self.tx_id].fusion:
            seq, variant_coordinates = self.get_variant_sequence_fusion(
                seq=tx_seq.seq, variants=variants
            )
        else:
            location = FeatureLocation(start=0, end=len(tx_seq.seq))
            seq, variant_coordinates = self.get_variant_sequence(
                seq=tx_seq.seq, location=location, offset=0,
                variants=variants[self.tx_id].transcriptional, pool=variants
            )

        stop_lost = self.check_stop_lost(seq, variant_coordinates)
        silent_mutation = self.check_silent_mutation(seq, variant_coordinates)

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
                last_m = aa_seq.seq[:lhs].rfind('M')
                if last_m > 0 and not tx_model.is_protein_coding:
                    actual_cds_start = cds_start + last_m * 3
                else:
                    actual_cds_start = cds_start

                prev_cds_start = self.find_prev_cds_start_same_frame(
                    cds_start=actual_cds_start,
                    cds_start_positions=cds_start_positions
                )

                for k in range(j + 1, min([j + 3, len(sites) - 1]) + 1):
                    rhs = sites[k]
                    tx_lhs = cds_start + lhs * 3
                    tx_rhs = cds_start + rhs * 3
                    if not self.has_any_variant(tx_lhs, tx_rhs, actual_cds_start,
                            prev_cds_start, variant_coordinates, stop_lost,
                            silent_mutation):
                        continue

                    peptides = [aa_seq[lhs:rhs]]
                    if peptides[0].seq.startswith('M') and lhs == 0:
                        peptides.append(peptides[0][1:])
                    for peptide in peptides:
                        is_valid = self.peptide_is_valid(peptide.seq)
                        if is_valid:
                            self.variant_peptides.add(str(peptide.seq))

    def generate_variant_comb(self) -> Iterable[seqvar.VariantRecordPool]:
        """ Generate combination of variants. """
        variant_type_mapper:Dict[seqvar.VariantRecord, str] = {}
        for variant in self.variant_pool[self.tx_id].transcriptional:
            variant_type_mapper[variant] = 'transcriptional'
        for variant in self.variant_pool[self.tx_id].intronic:
            variant_type_mapper[variant] = 'intronic'
        for variant in self.variant_pool[self.tx_id].fusion:
            variant_type_mapper[variant] = 'fusion'
            accepter_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
            if accepter_tx_id not in self.variant_pool:
                continue
            accepter_var_series = self.variant_pool[accepter_tx_id]
            for accepter_var in accepter_var_series.transcriptional:
                variant_type_mapper[accepter_var] = 'transcriptional'
            for accepter_var in accepter_var_series.intronic:
                variant_type_mapper[accepter_var] = 'intronic'
        for variant in self.variant_pool[self.tx_id].circ_rna:
            variant_type_mapper[variant] = 'circ_rna'

        all_variants = list(variant_type_mapper.keys())

        for i in range(len(all_variants)):
            for inds in combinations(range(len(all_variants)), i + 1):
                variants = [all_variants[i] for i in inds]
                pool = seqvar.VariantRecordPool()
                for variant in variants:
                    var_type = variant_type_mapper[variant]
                    tx_id = variant.transcript_id
                    if var_type == 'transcriptional':
                        pool.add_transcriptional_variant(variant, tx_id)
                    elif var_type == 'intronic':
                        pool.add_intronic_variant(variant, tx_id)
                    elif var_type == 'fusion':
                        pool.add_fusion_variant(variant, tx_id)
                    elif var_type == 'circ_rna':
                        pool.add_circ_rna(variant, tx_id)
                pool.sort()
                if self.has_incompatible_variants(pool):
                    continue
                yield pool


    def call_peptides(self):
        """ Call variant peptides """
        for comb in self.generate_variant_comb():
            if comb[self.tx_id].circ_rna:
                pass
            else:
                self.call_peptides_main(comb)

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
    variant_peptides:Set[str] = set()
    for tx_id in variant_pool.data.keys():
        caller = BruteForceVariantPeptideCaller()
        caller.variant_pool = variant_pool
        caller.reference_data = reference_data
        caller.tx_id = tx_id
        caller.tx_model = caller.reference_data.anno.transcripts[caller.tx_id]
        caller.tx_seq = caller.tx_model.get_transcript_sequence(
            caller.reference_data.genome['chr1']
        )

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

        caller.call_peptides()

        variant_peptides.update(caller.variant_peptides)

    for peptide in sorted(variant_peptides):
        print(peptide, file=sys.stdout)
