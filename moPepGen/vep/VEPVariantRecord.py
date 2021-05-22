""" This module defines data structure for variant records """
from __future__ import annotations
from typing import List
import re
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import dna
from moPepGen.vep.VepRecord import VEPRecord


_CONSEQUENCES_SILENT = [
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant'
]

class VariantRecord():
    """ Defines the location, ref and alt of a genomic variant.

    A VCF-like definition is used here to represent a variant. For any type of
    variant, both the ref and alt are at least length of one.

    Examples:
        * SNP/SNV:
          The below example represents a SNP/SNV at location 101, that a A is
          substituted with T.
            >>> VariantRecord(
                    location=FeatureLocation(100, 101),
                    ref=DNASeqRecord(Seq('A')),
                    alt=DNASeqRecord(Seq('T'))
                )
        
        * insertion
          The example below represents a insertion between 115 and 116, and the
          inserted sequence is CTGAA
            >>> VariantRecord(
                    location=FeatureLocation(115, 116),
                    ref='G',
                    alt='GCTGAA'
                )
        
        * deletion
          The example below represents a deletion of CAC from position 80 to
          82.
            >>> VariantRecord(
                    location=FeatureLocation(79, 83),
                    ref='CCAC',
                    alt='C'
                )
    
    Attributes:
        location (FeatureLocation): The location of the variant.
        ref (str): Reference sequence
        alt (str): Altered sequence
    """
    def __init__(self, location:FeatureLocation, ref:str, alt:str):
        """ Construct a VariantRecord object.
        
        Args:
            location (FeatureLocation): The location of the variant.
            ref (str): Reference sequence
            alt (str): Altered sequence
        """
        if len(location) != len(ref):
            raise ValueError("Length of ref must match with location.")
        self.location = location
        self.ref = ref
        self.alt = alt
    
    def __repr__(self) -> str:
        """Return representation of the VEP record."""
        return f"{self.location.start}:{self.location.end} {self.ref} ->" +\
            f" {self.alt}"


class VEPVariantRecord():
    """ A VEPVariantRecord object holds a variant record on a particular
    transcript. It contains the information of the variant including the
    transcript ID, and the location, ref, and alt of this transcript.
    
    Attributes:
        variant (VariatnRecord): the variant location at transcript
        consequence (List[str]): consequence type of this variant.
            See: https://uswest.ensembl.org/info/genome/variation/prediction/
            predicted_data.html#consequences
    
    Additional Attributes:
        chrom (str): The chromosome/scaffold seqname.
    """
    _silent_consequences = _CONSEQUENCES_SILENT
    def __init__(self, variant: VariantRecord, consequences: List[str]):
        """ Construct a VEPRecord object. """
        self.variant = variant
        self.consequences = consequences
    
    @property
    def transcript_id(self):
        """ Returns the transcirpt_id """
        return self.variant.location.seqname
    
    def __str__(self) -> str:
        """ return string """
        return f"{self.variant.location.start}-{self.variant.ref}" +\
            f":{self.variant.alt}"
    
    def __repr__(self) -> str:
        """Return representation of the VEP record."""
        consequences = '|'.join(self.consequences)
        return f"< {self.transcript_id}, {consequences}>"
    
    def __eq__(self, other: VEPVariantRecord)->bool:
        """ equal """
        return self.transcript_id == other.transcript_id and \
            self.variant.location == \
                other.variant.location and \
            self.variant.alt == other.variant.alt
        
    def __ne__(self, other: VEPVariantRecord) -> bool:
        """ not equal """
        return not self == other
    
    def __gt__(self, other: VEPVariantRecord) -> bool:
        """ greater than """
        if self.variant.location > other.variant.location:
            return True
        elif self.variant.location < other.variant.location:
            return False
        elif self.variant.alt > other.variant.alt:
            return True
        else:
            return False
    
    def __ge__(self, other: VEPVariantRecord) -> bool:
        """ greater or equal """
        return self > other or self == other
    
    def __lt__(self, other: VEPVariantRecord) -> bool:
        """ less than """
        return not (self > other and self == other)
    
    def __le__(self, other: VEPVariantRecord) -> bool:
        """ less or equal """
        return not (self > other)

    @staticmethod
    def from_vep_record(vep:VEPRecord, seq:dna.DNASeqRecord,
            transcript_id:str) -> VEPVariantRecord:
        """ It taks a VEPRecord and construct a TranscriptVariant Record object
        with only the variant information on the corresponding transcript.

        Args:
            vep (VEPRecord): The VEP record
            seq (DNASeqRecord): The transcript sequence
        """
        alt_position = vep.cdna_position.split('-')
        alt_start = int(alt_position[0]) - 1
        
        codon_ref, codon_alt = vep.codons

        if codon_ref == '-':
            alt_end = alt_start + 1
            ref = seq.seq[alt_start:alt_end]
            alt = ref + codon_alt
            alt_end = alt_start + len(ref)
        elif codon_alt == '-':
            alt_start -= 1
            alt_end = alt_start + len(codon_ref) + 1
            ref = seq.seq[alt_start:alt_end]
            alt = seq.seq[alt_start:alt_start+1]
        else:
            pattern = re.compile('[ATCG]+')
            match_ref = pattern.search(codon_ref)
            match_alt = pattern.search(codon_alt)
            if match_ref is not None:
                if match_alt is not None:
                    ref = match_ref.group()
                    alt_end = alt_start + len(ref)
                    alt = match_alt.group()
                else:
                    ref = match_ref.group()
                    alt_start -= 1
                    ref = seq.seq[alt_start] + ref
                    alt_end = alt_start + len(ref)
                    alt = seq.seq[alt_start:alt_end]
            elif match_alt is not None:
                alt = match_alt.group()
                # alt_start -= 1
                alt_end = alt_start + 1
                ref = seq.seq[alt_start:alt_end]
                alt = ref + alt
            else:
                raise ValueError('No alteration found in this VEP record')

        try:
            variant = VariantRecord(
                location=FeatureLocation(
                    seqname=transcript_id,
                    start=alt_start,
                    end=alt_end
                ),
                ref=ref,
                alt=alt
            )
        except ValueError as e:
            raise ValueError(e.args[0] + f' [{transcript_id}]')
        
        return VEPVariantRecord(
            variant=variant,
            consequences=vep.consequences
        )
    
    def is_silent(self) -> bool:
        """ Checks if consequence is silent """
        for consequence in self.consequences:
            if consequence in self._silent_consequences:
                return True
        return False
    
    def is_deletion(self) -> bool:
        """ Checks if is a deletion """
        return len(self.variant.ref) > 1 and len(self.variant.alt) == 1
    
    def is_insertion(self) -> bool:
        """ Checks if is a insertion """
        return len(self.variant.ref) == 1 and len(self.variant.alt) > 1

    def is_start_codon_retained(self,
            seq:dna.DNASeqRecordWithCoordinates) -> bool:
        """ Returns if the start codon is retained """
        start_codon = FeatureLocation(
            start=seq.orf.start,
            end=seq.orf.start + 3
        )
        # If the location of the variant does not overlaps with the start 
        # start codon, it is not considered.
        if not self.variant.location.overlaps(start_codon):
            return False
        if self.is_deletion():
            return False
        new_seq = seq.mutate(self)
        shift = len(self.variant.location.alt) - len(self.variant.location.ref)
        return new_seq[start_codon.start + shift:start_codon.end + shift] \
            == 'ATG'

    def is_start_codon_lost(self,
            seq:dna.DNASeqRecordWithCoordinates) -> bool:
        """ Returns if the start codon is lost """
        start_codon = FeatureLocation(
            start=seq.orf.start,
            end=seq.orf.start + 3
        )
        # If the location of the variant does not overlaps with the start 
        # start codon, it is not considered.
        if not self.variant.location.overlaps(start_codon):
            return False
        # If this is a deletion that the start codon is involved, treat it as
        # star codon lost, to search for new orf.
        if self.is_deletion():
            return True
        new_seq = seq.mutate(self)
        shift = len(self.variant.location.alt) - len(self.variant.location.ref)
        return new_seq[start_codon.start + shift:start_codon.end + shift] \
            != 'ATG'

    def is_frameshifting(self) -> bool:
        """ Retruns if is a frameshifting variant """
        return 'frameshift_variant' in self.consequences
    
    def is_cleave_site_gained(self, seq:DNASeqRecordWithCoordinates,
            rule:str, exception:str=None) -> bool:
        """ Returns if it gains a cleave site between the last and next cleave
        site. Only consider the case that a new cleave site is found between
        the last and next site. """
        left_codon_index = self.variant.location.start - \
            (self.variant.location.start - seq.orf.start) % 3
        right_codon_index = self.variant.location.end + 3 - \
            (self.variant.location.start + 3 - seq.orf.start) % 3
        left = seq.find_last_cleave_position(end=left_codon_index,
            rule=rule, exception=exception)
        right = seq.find_next_cleave_position(start=right_codon_index,
            rule=rule, exception=exception)
        return seq[left:right].mutate(self).iter_enzymatic_cleave_sites\
            is not None
    
    def is_cleave_site_lost(self, seq:dna.DNASeqRecordWithCoordinates,
            rule:str, exception:str=None) -> bool:
        """ Returns if it loses a cleave site """
        left_codon_index = self.variant.location.start - \
            (self.variant.location.start - seq.orf.start) % 3
        right_codon_index = self.variant.location.end + 3 - \
            (self.variant.location.start + 3 - seq.orf.start) % 3
        left_site = seq.find_cleave_positions_before(
            position=left_codon_index,
            rule=rule,
            exception=exception
        )
        right_site = seq.find_cleave_positions_after(
            position=right_codon_index,
            rule=rule,
            exception=exception
        )
        new_seq = seq.mutate(self)
        new_right_site = right_site + len(self.variant.alt) \
            - len(self.variant.ref)
        return seq[left_site:right_site].iter_enzymatic_cleave_sites() \
            is not None and \
            new_seq[left_site:new_right_site].iter_enzymatic_cleave_sites() \
            is None
        