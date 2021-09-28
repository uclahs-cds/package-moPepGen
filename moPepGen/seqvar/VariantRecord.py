""" Module for generic variant record """
from __future__ import annotations
import copy
from typing import TYPE_CHECKING
from moPepGen import ERROR_NO_TX_AVAILABLE, \
    ERROR_VARIANT_NOT_IN_GENE_COORDINATE, ERROR_INDEX_IN_INTRON, \
        ERROR_REF_LENGTH_NOT_MATCH_WITH_LOCATION
from moPepGen.SeqFeature import FeatureLocation


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.dna import DNASeqDict

_VARIANT_TYPES = ['SNV', 'INDEL', 'Fusion', 'RNAEditingSite',
    'Insertion', 'Deletion', 'Substitution', 'circRNA']
SINGLE_NUCLEOTIDE_SUBSTITUTION = ['SNV', 'SNP', 'INDEL', 'RNAEditingSite']
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
        if _type not in ['Substitution', 'Deletion', 'circRNA'] and \
                len(location) != len(ref):
            raise ValueError(ERROR_REF_LENGTH_NOT_MATCH_WITH_LOCATION)
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
        if self.type in ['Insertion', 'Substitution']:
            return int(self.attrs['DONOR_START'])
        raise ValueError(f"Don't know how to get donor start for variant type "
            f"{self.type}")

    def get_donor_end(self) -> int:
        """ Get donor end position """
        if self.type in ['Insertion', 'Substitution']:
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
        """ Convert to a GVF record. """
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

    def is_fusion(self) -> bool:
        """ Check if this is a fusion """
        return self.type == 'Fusion'

    def is_circ_rna(self) -> bool:
        """ check if this is a circRNA """
        return self.type == 'circRNA'

    def is_frameshifting(self) -> bool:
        """ Checks if the variant is frameshifting. """
        if self.type == 'Fusion':
            return True
        if self.type in ['Insertion', 'Substritution']:
            end = int(self.attrs['DONOR_END'])
            start = int(self.attrs['DONOR_START'])
            return abs(end - start) % 3 != 0
        if self.type == 'Deletion':
            return (self.location.end - self.location.start - 1) % 3 != 0
        return abs(len(self.alt) - len(self.ref)) % 3 != 0

    def frames_shifted(self):
        """ get number of nucleotide shifted """
        if self.type == 'Fusion':
            return 0
        if self.type in ['Insertion', 'Substritution']:
            end = int(self.attrs['DONOR_END'])
            start = int(self.attrs['DONOR_START'])
            return (end - start) % 3
        if self.type == 'Deletion':
            return (self.location.end - self.location.start - 1) % 3
        return (len(self.alt) - len(self.ref)) % 3

    def is_spanning_over_splicing_site(self, anno:GenomicAnnotation,
            transcript_id:str) -> bool:
        """ Check if this is spanning over splicing site """
        gene_id = self.location.seqname
        if gene_id not in anno.genes:
            raise ValueError(ERROR_VARIANT_NOT_IN_GENE_COORDINATE)
        start = self.location.start
        end = self.location.end

        try:
            anno.coordinate_gene_to_transcript(start, gene_id, transcript_id)
            start_in_intron = False
        except ValueError as e:
            if e.args[0] == ERROR_INDEX_IN_INTRON:
                start_in_intron = True
            else:
                raise e
        try:
            anno.coordinate_gene_to_transcript(end - 1, gene_id, transcript_id)
            end_in_intron = False
        except ValueError as e:
            if e.args[0] == ERROR_INDEX_IN_INTRON:
                end_in_intron = True
            else:
                raise e
        return ((not start_in_intron) and end_in_intron) \
                or (start_in_intron and (not end_in_intron))

    def to_transcript_variant(self, anno:GenomicAnnotation, genome:DNASeqDict,
            tx_id:str=None) -> VariantRecord:
        """ Get variant records with transcription coordinates """
        if not tx_id:
            tx_id = self.attrs['TRANSCRIPT_ID']
        if not tx_id:
            raise ValueError(ERROR_NO_TX_AVAILABLE)
        tx_model = anno.transcripts[tx_id]
        chrom = tx_model.transcript.chrom
        tx_seq = tx_model.get_transcript_sequence(genome[chrom])
        tx_start = tx_model.transcript.location.start
        start_genomic = anno.coordinate_gene_to_genomic(
            self.location.start, self.location.seqname
        )
        end_genomic = anno.coordinate_gene_to_genomic(
            self.location.end - 1, self.location.seqname
        )
        if tx_model.transcript.strand == -1:
            start_genomic, end_genomic = end_genomic, start_genomic
        end_genomic += 1
        if start_genomic < tx_start:
            if end_genomic < tx_start:
                raise ValueError(
                    'Variant not associated with the given transcript'
                )
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
                self.location.end - 1, self.location.seqname, tx_id
            ) + 1
            ref = self.ref
            alt = self.alt
            # If the distance from start to end on the transcript does not
            # equal to that on the gene, the variant contains at least one
            # intron. The ref is then regenerated.
            if end - start != len(self.location):
                ref = str(tx_seq.seq[start:end])
        location = FeatureLocation(seqname=tx_id, start=start, end=end)
        attrs = copy.deepcopy(self.attrs)
        del attrs['TRANSCRIPT_ID']
        attrs['GENE_ID'] = self.location.seqname
        return VariantRecord(
            location=location,
            ref=ref,
            alt=alt,
            _type=self.type,
            _id=self.id,
            attrs=attrs
        )

    def has_transcript(self):
        """ Checks if the variant has any transcript IDs annotated """
        return 'TRANSCRIPT_ID' in self.attrs
