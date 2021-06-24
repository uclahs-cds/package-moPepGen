""" Module for generic variant record """
from __future__ import annotations
from moPepGen.SeqFeature import FeatureLocation


_VARIANT_TYPES = ['SNV', 'INDEL', 'Fusion', 'RNAEditingSite',
    'AlternativeSplicing']
SINGLE_NUCLEOTIDE_SUBSTITUTION = ['SNV', 'SNP', 'INDEL']

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
    def __init__(self, location:FeatureLocation, ref:str, alt:str, _type:str,
            _id:str, attrs:dict=None):
        """ Construct a VariantRecord object.

        Args:
            location (FeatureLocation): The location of the variant.
            ref (str): Reference sequence
            alt (str): Altered sequence
            _type (str): Variant type, must be one of 'SNV', 'INDEL', 'Fusion',
                'RNAEditingSite', 'AlternativeSplicing'
        """
        if len(location) != len(ref):
            raise ValueError("Length of ref must match with location.")
        if _type not in _VARIANT_TYPES:
            raise ValueError('Variant type not supported')
        self.location = location
        self.ref = ref
        self.alt = alt
        self.type = _type
        self.id = _id
        self.attrs = attrs if attrs else {}

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

    def to_tvf(self) -> str:
        """ Convert to a TVF record. """
        chrom = self.location.seqname
        pos = str(int(self.location.start))
        _id = self.id
        qual = '.'
        _filter = '.'

        if self.type in SINGLE_NUCLEOTIDE_SUBSTITUTION:
            ref = str(self.ref)
            alt = str(self.alt)
        elif self.type == 'Fusion':
            ref = str(self.ref[0])
            alt = '<FUSION>'
        else:
            ref = str(self.ref[0])
            alt = f'<{self.type.upper()}>'

        info = self.info
        return '\t'.join([chrom, pos, _id, ref, alt, qual, _filter, info])

    @property
    def info(self) -> str:
        """ Get property of the INFO field """
        gene_id = self.attrs['GENE_ID']
        return f'GENE_ID={gene_id}'

    def is_snv(self) -> bool:
        """ Checks if the variant is a single nucleotide variant. """
        return len(self.ref) == 1 and len(self.alt) == 1

    def is_insertion(self) -> bool:
        """ Checks if the variant is an insertion """
        return len(self.ref) == 1 and len(self.alt) > 1

    def is_deletion(self) -> bool:
        """ Checks if the variant is a deletion. """
        return len(self.ref) > 1 and len(self.alt) == 1

    def is_frameshifting(self) -> bool:
        """ Checks if the variant is frameshifting. """
        return abs(len(self.alt) - len(self.ref)) % 3 != 0 or \
            self.type == 'Fusion'
