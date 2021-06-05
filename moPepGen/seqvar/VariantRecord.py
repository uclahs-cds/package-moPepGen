""""""
from __future__ import annotations
from moPepGen.SeqFeature import FeatureLocation


_VARIANT_TYPES = ['SNV', 'INDEL', 'Fusion', 'RNAEditingSite',
    'AlternativeSplicing']

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
    def __init__(self, location:FeatureLocation, ref:str, alt:str, type:str,
            id:str):
        f""" Construct a VariantRecord object.
        
        Args:
            location (FeatureLocation): The location of the variant.
            ref (str): Reference sequence
            alt (str): Altered sequence
            type (str): Variant type, must be one of ${_VARIANT_TYPES}
        """
        if len(location) != len(ref):
            raise ValueError("Length of ref must match with location.")
        if type not in _VARIANT_TYPES:
            raise ValueError('Variant type not supported')
        self.location = location
        self.ref = ref
        self.alt = alt
        self.type = type
        self.id = id
    
    def __hash__(self):
        """ hash """
        return hash((self.location.seqname, self.location.start,
             self.location.end, self.ref, self.alt, self.type))
    
    def __repr__(self) -> str:
        """Return representation of the VEP record."""
        return f"{self.location.start}:{self.location.end} {self.ref} ->" +\
            f" {self.alt}"
    
    def __eq__(self, other:VariantRecord) -> bool:
        """ equal to """
        return self.location == other.location
    
    def __ne__(self, other:VariantRecord) -> bool:
        """ not equal to """
        return not self == other
    
    def __gt__(self, other:VariantRecord) -> bool:
        """ greater than """
        return self.location > other.location
    
    def __ge__(self, other:VariantRecord) -> bool:
        """ greather or equal to """
        return self == other or self > other

    def __lt__(self, other:VariantRecord) -> bool:
        """ less then """
        return not self >= other
    
    def __le__(self, other:VariantRecord) -> bool:
        """ less or equal to """
        return not self > other

    def is_snv(self) -> bool:
        """"""
        return len(self.ref) == 1 and len(self.alt) == 1

    def is_insertion(self) -> bool:
        """"""
        return len(self.ref) == 1 and len(self.alt) > 1

    def is_deletion(self) -> bool:
        """"""
        return len(self.ref) > 1 and len(self.alt) == 1
    
    def is_frameshifting(self) -> bool:
        """"""
        return abs(len(self.alt) - len(self.ref)) % 3 != 0
