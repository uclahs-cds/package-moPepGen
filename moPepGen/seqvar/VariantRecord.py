""" Module for generic variant record """
from __future__ import annotations
import copy
from typing import List, TYPE_CHECKING
from moPepGen import ERROR_NO_TX_AVAILABLE
from moPepGen.SeqFeature import FeatureLocation


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.dna import DNASeqDict

_VARIANT_TYPES = ['SNV', 'INDEL', 'Fusion', 'RNAEditingSite',
    'Insertion', 'Deletion', 'Substitution', 'circRNA']
SINGLE_NUCLEOTIDE_SUBSTITUTION = ['SNV', 'SNP', 'INDEL']
ATTRS_POSITION = ['START', 'DONOR_START', 'ACCEPTER_START', 'ACCEPTER_POSITION']

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
                return self.type > other.type
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

    def get_donor_start(self) -> int:
        """ Get donor start position """
        if self.type == 'Insertion':
            return int(self.attrs['START'])
        if self.type == 'Substitution':
            return int(self.attrs['DONOR_START'])
        raise ValueError(f"Don't know how to get donor start for variant type "
            f"{self.type}")

    def get_donor_end(self) -> int:
        """ Get donor end position """
        if self.type == 'Insertion':
            return int(self.attrs['END'])
        if self.type == 'Substitution':
            return int(self.attrs['DONOR_END'])
        raise ValueError(f"Don't know how to get donor start for variant type "
            f"{self.type}")

    def get_accepter_position(self) -> int:
        """ Get accepter position for fusion only """
        if self.type != 'Fusion':
            raise ValueError("Don't know how to get accetoer start for "
                f"variant type {self.type}")
        return int(self.attrs['ACCEPTER_POSITION'])

    def get_alt_len(self) -> int:
        """ Get the length of the alt """
        if not self.alt.startswith('<'):
            return len(self.alt)
        if self.type == 'Fusion':
            return 1
        if self.type == 'Deletion':
            return 1
        if self.type in ['Insertion', 'Deletion']:
            return self.get_donor_end - self.get_donor_start()
        raise ValueError(f"Don't know how to get alt len for variant type "
            f"{self.type}")

    def to_string(self) -> str:
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
        elif self.type in ['Insertion', 'Deletion', 'Substitution']:
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
            if key in ATTRS_POSITION:
                val = str(int(val) + 1)
            elif isinstance(val, list):
                val = ','.join([str(x) for x in val])
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

    def to_transcript_variants(self, anno:GenomicAnnotation, genome:DNASeqDict,
            tx_ids:List[str]=None) -> List[VariantRecord]:
        """ Get variant records with transcription coordinates """
        if not tx_ids:
            tx_ids = self.attrs['TRANSCRIPTS']
        if not tx_ids:
            raise ValueError(ERROR_NO_TX_AVAILABLE)
        variants = []
        for tx_id in tx_ids:
            tx_model = anno.transcripts[tx_id]
            chrom = tx_model.transcript.chrom
            tx_seq = tx_model.get_transcript_sequence(genome[chrom])
            tx_start = tx_model.transcript.location.start
            start_genomic = anno.coordinate_gene_to_genomic(
                self.location.start, self.location.seqname
            )
            end_genomic = anno.coordinate_gene_to_genomic(
                self.location.end, self.location.seqname
            )
            if start_genomic < tx_start:
                if end_genomic < tx_start:
                    continue
                if not tx_model.is_cds_start_nf():
                    continue
                start = 0
                end = anno.coordinate_gene_to_transcript(
                    self.location.end + 1, self.location.seqname, tx_id
                )
                ref = str(tx_seq.seq[0:end])
                alt = str(self.alt[1:] + ref[-1])
            else:
                start = anno.coordinate_gene_to_transcript(
                    self.location.start, self.location.seqname, tx_id
                )
                end = anno.coordinate_gene_to_transcript(
                    self.location.end, self.location.seqname, tx_id
                )
                ref = self.ref
                alt = self.alt
            location = FeatureLocation(seqname=tx_id, start=start, end=end)
            attrs = copy.deepcopy(self.attrs)
            del attrs['TRANSCRIPTS']
            attrs['GENE_ID'] = self.location.seqname
            variant = VariantRecord(
                location=location,
                ref=ref,
                alt=alt,
                _type=self.type,
                _id=self.id,
                attrs=attrs
            )
            variants.append(variant)
        variants.sort()
        return variants

    def has_transcripts(self):
        """ Checks if the variant has any transcript IDs annotated """
        return 'TRANSCRIPTS' in self.attrs and len(self.attrs['TRANSCRIPTS']) > 0