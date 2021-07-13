""" Module for generic variant record """
from __future__ import annotations
from moPepGen.SeqFeature import FeatureLocation


_VARIANT_TYPES = ['SNV', 'INDEL', 'Fusion', 'RNAEditingSite',
    'Insertion', 'Deletion', 'Substitution']
SINGLE_NUCLEOTIDE_SUBSTITUTION = ['SNV', 'SNP', 'INDEL']
ATTRS_START = ['START', 'DONOR_START', 'ACCEPTOR_START']

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
        if _type not in ['Substitution', 'Deletion'] and \
                len(location) != len(ref):
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
        return self.location == other.location and \
            self.ref == other.ref and \
            self.type == other.type

    def __ne__(self, other:VariantRecord) -> bool:
        """ not equal to """
        return not self == other

    def __gt__(self, other:VariantRecord) -> bool:
        """ greater than """
        if self.location > other.location:
            return True
        if self.location == other.location:
            if self.alt > other.alt:
                return True
            if self.ref == other.ref:
                self.type > other.type
        return False

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
        # using 1-base position
        pos = str(int(self.location.start) + 1)
        _id = self.id
        qual = '.'
        _filter = '.'

        if self.type in SINGLE_NUCLEOTIDE_SUBSTITUTION:
            ref = str(self.ref)
            alt = str(self.alt)
        elif self.type == 'Fusion':
            ref = str(self.ref[0])
            alt = '<FUSION>'
        elif self.type in ATTRS_START:
            ref = str(self.ref[0])
            alt = f'<{self.type.upper()[:3]}>'
        else:
            ref = str(self.ref[0])
            alt = f'<{self.type.upper()}>'

        info = self.info
        return '\t'.join([chrom, pos, _id, ref, alt, qual, _filter, info])

    @property
    def info(self) -> str:
        """ Get property of the INFO field """
        out = ''
        for key,val in self.attrs.items():
            # using 1-base position
            if key in ['DONOR_START', 'ACCEPTOR_START', 'START']:
                val = str(int(val) + 1)
            out += f'{key.upper()}={val};'
        return out.rstrip(';')

    def is_snv(self) -> bool:
        """ Checks if the variant is a single nucleotide variant. """
        return len(self.ref) == 1 and len(self.alt) == 1

    def is_insertion(self) -> bool:
        """ Checks if the variant is an insertion """
        if self.type == 'Insertion':
            return True
        return len(self.ref) == 1 and len(self.alt) > 1

    def is_deletion(self) -> bool:
        """ Checks if the variant is a deletion. """
        if self.type == 'Deletion':
            return True
        return len(self.ref) > 1 and len(self.alt) == 1

    def is_frameshifting(self) -> bool:
        """ Checks if the variant is frameshifting. """
        if self.type == 'Fusion':
            return True
        if self.type == 'Insertion':
            end = int(self.attrs['END'])
            start = int(self.attrs['START'])
            return abs(end - start) % 3 != 0
        if self.type == 'Substritution':
            end = int(self.attrs['DONOR_END'])
            start = int(self.attrs['DONOR_START'])
            return abs(end - start) % 3 != 0
        if self.type == 'Deletion':
            return (self.location.end - self.location.start - 1) % 3 != 0
        return abs(len(self.alt) - len(self.ref)) % 3 != 0
