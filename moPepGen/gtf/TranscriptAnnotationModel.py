""" Module for Transcript Annotation model """
from __future__ import annotations
from typing import List, TYPE_CHECKING
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen import dna, ERROR_INDEX_IN_INTRON


if TYPE_CHECKING:
    from Bio.Seq import Seq

GTF_FEATURE_TYPES = {
    'transcript': 'transcript',
    'cds': 'cds',
    'exon': 'exon',
    'start_codon': 'start_codon',
    'stop_codon': 'stop_codon',
    'utr': 'utr',
    'selenocysteine': 'selenocysteine',
    'five_prime_utr': 'five_utr',
    'three_prime_utr': 'three_utr'
}

class TranscriptAnnotationModel():
    """ A TranscriptAnnotationModel holds all the annotations associated with
    the same transcript from a GTF file.

    Attributes:
        transcript (GTFRecord): The GTF record of the transcript
        cds (List[GTFRecord]): The GTF records of all the coding sequences of
            the transcript.
        exon (List[GTFRecord]): The GTF records of all the exons of the
            transcript.
        start_codon (List[GTFRecord]): The GTF records of all the start codons
            of the transcript.
        stop_codon (List[GTFRecord]): The GTF records of all the stop codons
            of the transcript.
    """
    def __init__(self, transcript:GTFSeqFeature=None,
            cds:List[GTFSeqFeature]=None, exon:List[GTFSeqFeature]=None,
            start_codon:List[GTFSeqFeature]=None,
            stop_codon:List[GTFSeqFeature]=None,
            utr:List[GTFSeqFeature]=None,
            five_utr:List[GTFSeqFeature]=None,
            three_utr:List[GTFSeqFeature]=None,
            selenocysteine:List[GTFSeqFeature]=None,
            is_protein_coding:bool=None,
            _seq:dna.DNASeqRecordWithCoordinates=None,
            transcript_id:str=None,
            gene_id:str=None,
            protein_id:str=None,
            gene_name:str=None,
            gene_type:str=None):
        """ Construct a TranscriptAnnotationmodel """
        self.transcript = transcript
        self.cds = cds or []
        self.exon = exon or []
        self.start_codon = start_codon or []
        self.stop_codon = stop_codon or []
        self.utr = utr or []
        self.five_utr = five_utr or []
        self.three_utr = three_utr or []
        self.selenocysteine = selenocysteine or []
        self.is_protein_coding = is_protein_coding
        self._seq = _seq
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.protein_id = protein_id
        self.gene_name = gene_name
        self.gene_type = gene_type

    def add_record(self, _type:str, record: GTFSeqFeature):
        """ Add a GTFRecrod into a TranscriptAnnotationModel. If the biotype
        is cds, exon, start_codon, or stop_codon, it is inserted in order.

        Args:
            type (str): Type of annotation to add. Must be from transcript,
                cds, exon, start_codon, or stop_codon.
            record (GTFRecord): The GTF record to be added.
        """
        if _type not in GTF_FEATURE_TYPES:
            raise ValueError(f'Type must be from {list(GTF_FEATURE_TYPES.keys())}')
        attr = GTF_FEATURE_TYPES[_type]

        for key in ['transcript_id', 'gene_id', 'protein_id', 'gene_name', 'gene_type']:
            if hasattr(record, key):
                if self.__getattribute__(key) is None:
                    val = record.__getattribute__(key)
                    if val is not None and val != '':
                        self.__setattr__(key, record.__getattribute__(key))
                else:
                    record.__setattr__(key, self.__getattribute__(key))

        if _type == 'transcript':
            self.transcript = record
            if 'is_protein_coding' in record.attributes:
                is_protein_coding = record.attributes.pop('is_protein_coding')
                self.is_protein_coding = is_protein_coding == 'true'
        else:
            if self.__getattribute__(attr) is None:
                self.__setattr__(attr, [])
            self.__getattribute__(attr).append(record)

    def split_utr(self) -> None:
        """ Infer 3'UTR or 5'UTR. CDS must be already sorted. """
        if not self.utr:
            return
        if not self.cds:
            raise ValueError(
                "UTR found but not CDS for transcript" + \
                    self.transcript.transcript_id
            )
        strand = self.transcript.strand
        for utr in self.utr:
            if strand == 1 and utr < self.cds[0] \
                    or strand == -1 and utr > self.cds[-1]:
                self.five_utr.append(utr)
            else:
                self.three_utr.append(utr)

    def remove_cached_seq(self):
        """ Remove the cased sequence """
        self._seq = None

    def sort_records(self):
        """ sort records """
        self.cds.sort()
        self.exon.sort()
        self.start_codon.sort()
        self.stop_codon.sort()
        self.split_utr()
        self.utr.sort()
        self.three_utr.sort()
        self.five_utr.sort()
        self.selenocysteine.sort()

    def get_cds_start_index(self) -> int:
        """ Returns the CDS start index of the transcript """
        cds_start = 0
        if self.transcript.strand == 1:
            cds_frame = self.cds[0].frame
            for exon in self.exon:
                if self.cds[0].location.start not in exon:
                    cds_start += len(exon)
                    continue
                cds_start += (self.cds[0].location.start - exon.location.start)
                break
        elif self.transcript.strand == -1:
            cds_frame = self.cds[-1].frame or 0
            for exon in reversed(self.exon):
                if self.cds[-1].location.end != exon.location.end and \
                        self.cds[-1].location.end not in exon:
                    cds_start += len(exon)
                    continue
                cds_start += (exon.location.end - self.cds[-1].location.end)
                break
        else:
            raise ValueError('Strand must not be unknown.')
        return cds_start + cds_frame

    def get_cds_end_index(self, seq:Seq, start:int) -> int:
        """ Returns the CDS stop index of the transcript. """
        if self.three_utr:
            if self.transcript.strand == 1:
                end = self.get_transcript_index(self.three_utr[0].location.start)
            else:
                end = self.get_transcript_index(self.three_utr[-1].location.end - 1)
            return end - (end - start) % 3
        return len(seq) - (len(seq) - start) % 3

    def get_transcript_sequence(self, chrom:dna.DNASeqRecord,
            cache:bool=False) -> dna.DNASeqRecordWithCoordinates:
        """ Returns the transcript sequence. The is done by concating all the
        exon sequences. If the gene is on the negatice strand, the reverse
        complement is returned.

        Args:
            chrom (DNASeqRecord): The chromosome sequence that the transcript
                is located.
        """
        if self._seq:
            return self._seq
        if len(self.exon) == 0:
            raise ValueError("Transcript model has no exon")
        seq = None
        for location in [it.location for it in self.exon]:
            new_seq = chrom.seq[location.start:location.end]
            if seq is None:
                seq = new_seq
            else:
                seq = seq + new_seq

        if self.transcript.strand == -1:
            seq = seq.reverse_complement()

        if self.cds:
            cds_start = self.get_cds_start_index()
            cds_end = self.get_cds_end_index(seq, cds_start)
            orf = FeatureLocation(start=cds_start, end=cds_end)
        else:
            orf = None

        selenocystein = []
        for sec in self.selenocysteine:
            if sec.strand == 1:
                sec_start = self.get_transcript_index(sec.location.start)
                sec_end = self.get_transcript_index(sec.location.end - 1) + 1
            else:
                sec_start = self.get_transcript_index(sec.location.end - 1)
                sec_end = self.get_transcript_index(sec.location.start) + 1
            selenocystein.append(FeatureLocation(start=sec_start, end=sec_end))
        selenocystein.sort()

        location = MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.transcript_id,
                start=0,
                end=len(seq)
            )
        )
        transcript = dna.DNASeqRecordWithCoordinates(
            seq=seq,
            orf=orf,
            selenocysteine=selenocystein,
            locations=[location]
        )
        transcript.id = self.transcript.transcript_id
        transcript.name = self.transcript.transcript_id
        transcript.description = self.transcript.transcript_id +\
            '|' + self.transcript.gene_id
        if 'protein_id' in self.transcript.attributes:
            transcript.description += '|' + \
                self.transcript.protein_id
        if cache:
            self._seq = transcript
        return transcript

    def get_cdna_sequence(self, chrom:dna.DNASeqRecord
            ) -> dna.DNASeqRecordWithCoordinates:
        """ Get the cDNA sequence.

        Args:
            chrom (DNASeqRecord): The chromosome sequence.

        Returns:
            The cDNA sequence (DNASeqRecordWithCoordinates)
        """
        if len(self.cds) == 0:
            raise ValueError('Transcript model has no cds')
        seq = None
        locations = [it.location for it in self.cds]
        for location in locations:
            new_seq = chrom.seq[location.start:location.end]
            if seq is None:
                seq = new_seq
            else:
                seq = seq + new_seq

        cds_start = self.get_cds_start_index()

        if self.transcript.strand == -1:
            seq = seq.reverse_complement()

        location = MatchedLocation(
            query=FeatureLocation(start=0, end=len(seq)),
            ref=FeatureLocation(
                seqname=self.transcript.transcript_id,
                start=cds_start,
                end=cds_start + len(seq)
            )
        )
        cdna = dna.DNASeqRecordWithCoordinates(
            seq=seq,
            locations=[location]
        )
        cdna.id = self.transcript.transcript_id
        cdna.name = self.transcript.transcript_id
        cdna.description = self.transcript.transcript_id + '|' \
            + self.transcript.gene_id + '|' \
            + self.transcript.protein_id
        return cdna

    def get_transcript_index(self, genomic_index:int) -> int:
        """ Get the transcript index from a genomic index. """
        if genomic_index < self.exon[0].location.start \
                or genomic_index >= self.exon[-1].location.end:
            raise ValueError(
                r'The genomic index isn\'t in the range of this transcript'
            )
        if self.transcript.strand == 1:
            index = 0
            for exon in self.exon:
                if exon.location.end < genomic_index:
                    index += exon.location.end - exon.location.start
                elif exon.location.end == genomic_index:
                    raise ValueError(ERROR_INDEX_IN_INTRON)
                elif exon.location.start <= genomic_index:
                    index += genomic_index - exon.location.start
                    break
                else:
                    raise ValueError(ERROR_INDEX_IN_INTRON)
        elif self.transcript.strand == -1:
            index = -1
            for exon in reversed(self.exon):
                if exon.location.start >= genomic_index:
                    index += exon.location.end - exon.location.start
                    if exon.location.start == genomic_index:
                        break
                elif exon.location.end > genomic_index:
                    index += exon.location.end - genomic_index
                    break
                else:
                    raise ValueError(ERROR_INDEX_IN_INTRON)
        return index

    def get_transcript_start_genomic_coordinate(self) -> int:
        """ Get the genomic coordinate of the start point of the transcript """
        if self.transcript.location.strand == 1:
            return int(self.exon[0].location.start)
        return int(self.exon[-1].location.end)

    def is_cds_start_nf(self) -> bool:
        """ Returns if the transcript has the tag of cds_start_NF """
        return 'tag' in self.transcript.attributes and \
            'cds_start_NF' in self.transcript.attributes['tag']

    def is_mrna_end_nf(self) -> bool:
        """ Returns if the transcript has the tag of mrna_end_NF """
        return 'tag' in self.transcript.attributes and \
            'mRNA_end_NF' in self.transcript.attributes['tag']

    def transcript_len(self) -> int:
        """ Get the transcript length minus introns """
        length = 0
        for exon in self.exon:
            length += exon.location.end - exon.location.start
        return length

    def is_exonic(self, pos:int) -> bool:
        """ Checks if the given genomic position is contained by any exon of
        the current transcript model.

        Args:
            pos (int): A genomic position
        """
        for exon in self.exon:
            if exon.location.start <= pos < exon.location.end:
                return True
        return False

    def get_next_exon(self, pos:int) -> GTFSeqFeature:
        """ Get the next exon downstream to the position. The position given
        must be genomic position """
        if self.transcript.strand == 1:
            for exon in self.exon:
                if exon.location.start > pos:
                    return exon
        else:
            for exon in reversed(self.transcript.exon):
                if exon.location.start <= pos:
                    return exon
        return None

    def get_upstream_exon_end(self, pos:int) -> int:
        """ Find the upstream exon end """

        if self.transcript.strand == 1:
            for exon in self.exon:
                if exon.location.end > pos:
                    break
                ind = exon.location.end - 1
        else:
            for exon in reversed(self.exon):
                if exon.location.start < pos:
                    break
                ind = exon.location.start
        if ind == -1:
            raise ValueError("Could not find the upstream exon end.")
        return ind

    def get_downstream_exon_start(self, pos:int) -> int:
        """ Find the downstream exon end """
        ind = -1
        if self.transcript.strand == 1:
            for exon in self.exon:
                if exon.location.start >= pos:
                    ind = exon.location.start
                    break
        else:
            for exon in reversed(self.exon):
                if exon.location.end - 1 <= pos:
                    ind = exon.location.end - 1
                    break
        if ind == -1:
            raise ValueError("Could not find the upstream exon end.")
        return ind

    def get_exon_with_start(self, start:int, offset:int=0) -> int:
        """ Get the exon with a given START position. """
        for i, exon in enumerate(self.exon[offset:]):
            if exon.location.start == start:
                return i + offset
        return -1

    def get_exon_with_end(self, end:int, offset:int=0) -> int:
        """ Get the exon with a given end position. """
        for i, exon in enumerate(self.exon[offset:]):
            if exon.location.end == end:
                return i + offset
        return -1

    def get_exon_overlaps(self, loc:FeatureLocation, offset:int=0) -> int:
        """ Get exons that overlap with a given location """
        indices = []
        for i, exon in enumerate(self.exon[offset:]):
            if exon.location.overlaps(loc):
                indices.append(i + offset)
        return indices

    def get_exon_containing(self, pos:int) -> int:
        """ Get exon that contains a position """
        for i, exon in enumerate(self.exon):
            if pos in exon.location:
                return i
        return -1

    def get_exon_inner(self, loc:FeatureLocation, offset:int=0) -> int:
        """ get exon that are inner subset of a given location """
        indices = []
        for i, exon in enumerate(self.exon[offset:]):
            if loc.is_superset(exon.location):
                indices.append(i + offset)
        return indices

    def has_junction(self, junction:FeatureLocation) -> bool:
        """ Checks if has the given junction """
        for i, exon1 in enumerate(self.exon):
            if exon1 is self.exon[-1]:
                break
            exon2 = self.exon[i + 1]
            if exon2.location.start > junction.end:
                break
            if exon1.location.end == junction.start \
                    and exon2.location.start == junction.end:
                return True
        return False
