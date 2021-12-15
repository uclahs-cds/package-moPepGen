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
    from moPepGen.dna import DNASeqDict, DNASeqRecord

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

    @property
    def transcript_id(self) -> str:
        """ Transcript ID """
        if 'TRANSCRIPT_ID' in self.attrs:
            return self.attrs['TRANSCRIPT_ID']
        return self.location.seqname

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
        ref_len = len(self.location)
        if self.type in ['Insertion', 'Substritution']:
            end = int(self.attrs['DONOR_END'])
            start = int(self.attrs['DONOR_START'])
            alt_len = end - start
        elif self.type == 'Deletion':
            alt_len = 1
        else:
            alt_len = len(self.alt)
        return (ref_len - alt_len) % 3 != 0

    def frames_shifted(self):
        """ get number of nucleotide shifted
        TODO: tests needed """
        if self.type == 'Fusion':
            return True
        ref_len = len(self.location)
        if self.type in ['Insertion', 'Substritution']:
            end = int(self.attrs['DONOR_END'])
            start = int(self.attrs['DONOR_START'])
            alt_len = end - start
        elif self.type == 'Deletion':
            alt_len = 1
        else:
            alt_len = len(self.alt)
        return (ref_len - alt_len) % 3

    def is_spanning_over_splicing_site(self, anno:GenomicAnnotation,
            transcript_id:str) -> bool:
        """ Check if this is spanning over splicing site """
        gene_id = self.location.seqname
        if gene_id not in anno.genes:
            raise ValueError(ERROR_VARIANT_NOT_IN_GENE_COORDINATE)
        if self.is_fusion():
            start = self.location.start - 1
            end = self.location.end - 1
        else:
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
        gene_id = self.location.seqname
        strand = tx_model.transcript.strand

        tx_start = tx_model.transcript.location.start if strand == 1 else \
            tx_model.transcript.location.end - 1

        var_start = self.location.start
        var_end = self.location.end

        start_genomic = anno.coordinate_gene_to_genomic(var_start, gene_id)
        end_genomic = anno.coordinate_gene_to_genomic(var_end - 1, gene_id)
        end_genomic += 1

        variant_start_before_tx_start = \
            (strand == 1 and start_genomic < tx_start) or \
            (strand == -1 and tx_start < start_genomic)
        variant_end_before_tx_start = \
            (strand == 1 and end_genomic <= tx_start) or \
            (strand == -1 and tx_start <= end_genomic)

        if variant_start_before_tx_start:
            if variant_end_before_tx_start:
                raise ValueError(
                    'Variant not associated with the given transcript'
                )
            start = 0
            end = anno.coordinate_gene_to_transcript(var_end, gene_id, tx_id) + 1
            ref = str(tx_seq.seq[0:end])
            alt = str(self.alt[1:] + ref[-1])
        else:
            if self.is_fusion():
                start = anno.coordinate_gene_to_transcript(var_start - 1, gene_id, tx_id) + 1
                end = start + 1
                # if start == len(tx_seq):
                #     raise err.FusionBreakpointIsEndOfTranscript(self.id)
            else:
                start = anno.coordinate_gene_to_transcript(var_start, gene_id, tx_id)
                end = anno.coordinate_gene_to_transcript(var_end - 1, gene_id, tx_id) + 1
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

    def is_stop_lost(self, stop:int):
        """ Check if this is a stop lost mutation """
        loc = FeatureLocation(start=stop, end=stop+3)
        return self.location.overlaps(loc)

    def to_end_inclusion(self, seq:DNASeqRecord):
        """ Convert the variant to start exlusion and end inclusion format """
        location = FeatureLocation(
            seqname=self.location.seqname,
            start=self.location.start + 1,
            end=self.location.end + 1
        )
        ref = self.ref[1:] + str(seq.seq[self.location.end])
        alt = self.alt[1:] + str(seq.seq[self.location.end])
        self.location = location
        self.ref = ref
        self.alt = alt

    def shift_breakpoint_to_closest_exon(self, anno:GenomicAnnotation):
        """ Shift fusion breakpoints to the closest exon. Donor breakpoint will
        shift to the upstream and accepter breakpoint to the downstream. """
        if not self.is_fusion():
            raise ValueError(
                "Don't know how to shift breakpoint for non-fusion variants."
            )
        donor_gene_id = self.location.seqname
        donor_tx_id = self.attrs['TRANSCRIPT_ID']
        donor_tx_model = anno.transcripts[donor_tx_id]
        left_breakpoint = anno.coordinate_gene_to_genomic(
            index=self.location.start - 1, gene=donor_gene_id
        )
        if donor_tx_model.is_exonic(left_breakpoint):
            left_insertion_start = None
            left_insertion_end = None
        else:
            upstream_exon_end = donor_tx_model.get_upstream_exon_end(left_breakpoint)
            left_insertion_start = anno.coordinate_genomic_to_gene(
                index=upstream_exon_end, gene=donor_gene_id
            ) + 1
            left_insertion_end = self.location.start
            self.location = FeatureLocation(
                start=left_insertion_start,
                end=left_insertion_start + 1,
                seqname = donor_gene_id
            )

        accepter_gene_id = self.attrs['ACCEPTER_GENE_ID']
        accepter_tx_id = self.attrs['ACCEPTER_TRANSCRIPT_ID']
        accepter_tx_model = anno.transcripts[accepter_tx_id]
        right_breakpoint = anno.coordinate_gene_to_genomic(
            index=self.get_accepter_position(), gene=accepter_gene_id
        )
        if accepter_tx_model.is_exonic(right_breakpoint):
            right_insertion_start = None
            right_insertion_end = None
        else:
            downstream_exon_start = accepter_tx_model.get_downstream_exon_start(
                pos=right_breakpoint
            )
            right_insertion_start = self.get_accepter_position()
            right_insertion_end = anno.coordinate_genomic_to_gene(
                downstream_exon_start, accepter_gene_id
            )
            self.attrs['ACCEPTER_POSITION'] = anno.coordinate_genomic_to_gene(
                index=downstream_exon_start, gene=accepter_gene_id
            )

        self.attrs.update({
            'LEFT_INSERTION_START': left_insertion_start,
            'LEFT_INSERTION_END': left_insertion_end,
            'RIGHT_INSERTION_START': right_insertion_start,
            'RIGHT_INSERTION_END': right_insertion_end
        })
