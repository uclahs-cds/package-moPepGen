""" Module for Arriba parser """
from __future__ import annotations
import itertools
from typing import Iterable, List, IO
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, gtf, dna, err


def parse(handle:IO) -> Iterable[ArribaRecord]:
    """ parse Arriba's output """
    for line in handle:
        line:str
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        fields[9] = int(fields[9])
        fields[10] = int(fields[10])
        fields[11] = int(fields[11])
        fields[12] = int(fields[12])
        fields[13] = int(fields[13])
        fields[14] = ArribaConfidence(fields[14])
        fields[26] = fields[26].split(',')
        fields[29] = fields[29].split(',')

        yield ArribaRecord(*fields)

class ArribaConfidence():
    """ Arriba's confidence """
    levels = {'low': 0, 'medium': 1, 'high': 2}
    """ Arriba's confidence """
    def __init__(self, data:str):
        """ constructor """
        if data not in self.levels:
            raise ValueError(f"Invalid confidence of {data}")
        self.data = data

    def to_int(self) -> int:
        """ Convert to int """
        return self.levels[self.data]

    def __lt__(self, other:ArribaConfidence) -> bool:
        """ larger than """
        return self.to_int() > other.to_int()

    def __le__(self, other:ArribaConfidence) -> bool:
        """ larger or equal to """
        return self == other or self > other

    def __st__(self, other:ArribaConfidence) -> bool:
        """ smaller than """
        return not self >= other

    def __se__(self, other:ArribaConfidence) -> bool:
        """ smaller or equal to """
        return not self > other

    def __eq__(self, other:ArribaConfidence) -> bool:
        """ equal to """
        return self.data == other.data

    def __ne__(self, other:ArribaConfidence) -> bool:
        """ not equal to """
        return not self == other

class ArribaRecord():
    """ Arriba Record. More info see:
    https://arriba.readthedocs.io/en/latest/output-files/ """
    def __init__(self, gene1:str, gene2:str, strand1:str, strand2:str,
            breakpoint1:str, breakpoint2:str, site1:str, site2:str, _type:str,
            split_reads1:str, split_reads2:str, discordant_mates:int,
            coverage1:int, coverage2:int, confidence:ArribaConfidence,
            reading_frame:str, tags:str, retained_protein_domains:str,
            closest_genomic_breakpoint1:str, closest_genomic_breakpoint2:str,
            gene_id1:str, gene_id2:str, transcript_id1:str, transcript_id2:str,
            direction1:str, direction2:str, filters:List[str],
            fusion_transcript:str, peptide_sequence:str,
            read_identifiers:List[str]):
        """ Constructor """
        self.gene1 = gene1
        self.gene2 = gene2
        self.strand1 = strand1
        self.strand2 = strand2
        self.breakpoint1 = breakpoint1
        self.breakpoint2 = breakpoint2
        self.site1 = site1
        self.site2 = site2
        self._type = _type
        self.split_reads1 = split_reads1
        self.split_reads2 = split_reads2
        self.discordant_mates = discordant_mates
        self.coverage1 = coverage1
        self.coverage2 = coverage2
        self.confidence = confidence
        self.reading_frame = reading_frame
        self.tags = tags
        self.retained_protein_domains = retained_protein_domains
        self.closest_genomic_breakpoint1 = closest_genomic_breakpoint1
        self.closest_genomic_breakpoint2 = closest_genomic_breakpoint2
        self.gene_id1 = gene_id1
        self.gene_id2 = gene_id2
        self.transcript_id1 = transcript_id1
        self.transcript_id2 = transcript_id2
        self.direction1 = direction1
        self.direction2 = direction2
        self.filters = filters
        self.fusion_transcript = fusion_transcript
        self.peptide_sequence = peptide_sequence
        self.read_identifiers = read_identifiers

    @property
    def gene_strand1(self) -> int:
        """ gene strand1 """
        return self.infer_strand(self.strand1.split('/')[0])

    @property
    def gene_strand2(self) -> int:
        """ gene strand2 """
        return self.infer_strand(self.strand2.split('/')[0])

    @property
    def transcript_strand1(self) -> int:
        """ transcript strand1 """
        return self.infer_strand(self.strand1.split('/')[1])

    @property
    def transcript_strand2(self) -> int:
        """ transcript strand2 """
        return self.infer_strand(self.strand2.split('/')[1])

    @property
    def breakpoint1_position(self) -> int:
        """ Get breakpoint1 position """
        return int(self.breakpoint1.split(':')[1])

    @property
    def breakpoint2_position(self) -> int:
        """ Get breakpoint2 position """
        return int(self.breakpoint2.split(':')[1])

    @staticmethod
    def infer_strand(strand:str):
        """ Infer strand """
        if strand == '+':
            return 1
        if strand == '-':
            return -1
        if strand == '.':
            return 0
        raise ValueError(f"Can not infer strand from {strand}")

    def transcript_on_antisense_strand(self, anno:gtf.GenomicAnnotation) -> bool:
        """ Checks whether donor or accepter transcripts are on the antisense
        strand (-) or if the strand can not be interpreted (.) by Arriba. """
        gene_model1 = anno.genes[self.gene_id1]
        gene_model2 = anno.genes[self.gene_id2]
        return self.transcript_strand1 != gene_model1.strand or\
            self.transcript_strand2 != gene_model2.strand

    def get_donor_transcripts(self, anno:gtf.GenomicAnnotation
            ) -> List[gtf.TranscriptAnnotationModel]:
        """ Get all possible donor transcripts """
        pos = self.breakpoint1_position - 1
        return anno.get_transcripts_with_position(self.gene_id1, pos)

    def get_accepter_transcripts(self, anno:gtf.GenomicAnnotation
            ) -> List[gtf.TranscriptAnnotationModel]:
        """ Get all possible donor transcripts """
        pos = self.breakpoint2_position - 1
        return anno.get_transcripts_with_position(self.gene_id2, pos)

    def is_valid(self, min_split_reads1:int, min_split_reads2, confidence:str
            ) -> bool:
        """ Checks if the ArribaRecord is valid """
        return self.split_reads1 >= min_split_reads1 and \
            self.split_reads2 >= min_split_reads2 and \
            self.confidence >= ArribaConfidence(confidence)

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> List[seqvar.VariantRecord]:
        """ Convert the Arriba fusion record to a VariantRecord.
        Arriba uses a 1-based coordinate system.

        For Arriba, `breakpoint1` is used to indicate the genomic location of
        the breakpoint of the donor gene, while `breakpoint2` is to that of the
        accepter gene. `direction1` and `direction2` are either 'upstream' or
        'downstream' depends on donor/accepter and its strand. For example,
        for a donor gene on + strand, the value of `direction1` would be
        downstream meaning that the breakpoint is to the downstream side of
        `breakpoint1`. While an acceptor gene on - strand will have the
        `direction2` value being `downstream`. Noted that on genetic
        coordinates, the breakpoint of donor gene will always be 'downstream'
        and that of an accepter gene will always be 'upstream'.

        Args:
            anno (gtf.GenomicAnnotation): genomic annotation loaded from a GTF
                file.
            genome (dna.DNASeqDict): genome sequences loaded from a genome
                FASTA file.

        Return:
            A list of seqvar.VaraintRecord.
        """
        try:
            donor_gene_model = anno.genes[self.gene_id1]
        except KeyError as error:
            raise err.GeneNotFoundError(self.gene_id1) from error
        try:
            accepter_gene_model = anno.genes[self.gene_id2]
        except KeyError as error:
            raise err.GeneNotFoundError(self.gene_id2) from error

        donor_gene_symbol = donor_gene_model.gene_name
        donor_chrom = donor_gene_model.chrom
        left_breakpoint = self.breakpoint1_position
        donor_genome_position = f'{donor_chrom}:{left_breakpoint}:{left_breakpoint}'
        donor_position = anno.coordinate_genomic_to_gene(left_breakpoint - 1, self.gene_id1) + 1
        donor_transcripts = self.get_donor_transcripts(anno)

        accepter_gene_symbol = accepter_gene_model.gene_name
        accepter_chrom = accepter_gene_model.chrom
        right_breakpoint = self.breakpoint2_position
        accepter_genome_position = f'{accepter_chrom}:{right_breakpoint}:{right_breakpoint}'
        accepter_position = anno.coordinate_genomic_to_gene(right_breakpoint - 1, self.gene_id2)
        accepter_transcripts = self.get_accepter_transcripts(anno)

        records = []

        if donor_gene_model.strand == 1:
            ref_seq = genome[donor_chrom].seq[left_breakpoint]
        else:
            ref_seq = genome[donor_chrom]\
                .seq[left_breakpoint:left_breakpoint + 1]\
                .reverse_complement()
            ref_seq = str(ref_seq)

        perms = itertools.product(donor_transcripts, accepter_transcripts)
        for donor_tx, accepter_tx in perms:
            donor_tx_id = donor_tx.transcript.transcript_id
            accepter_tx_id = accepter_tx.transcript.transcript_id

            fusion_id = f'FUSION-{donor_tx_id}:{donor_position}'\
                f'-{accepter_tx_id}:{accepter_position}'

            location = FeatureLocation(
                seqname=self.gene_id1,
                start=donor_position,
                end=donor_position + 1
            )
            attrs = {
                'TRANSCRIPT_ID': donor_tx_id,
                'GENE_SYMBOL': donor_gene_symbol,
                'GENOMIC_POSITION': donor_genome_position,
                'ACCEPTER_GENE_ID': self.gene_id2,
                'ACCEPTER_TRANSCRIPT_ID': accepter_tx_id,
                'ACCEPTER_SYMBOL': accepter_gene_symbol,
                'ACCEPTER_POSITION': accepter_position,
                'ACCEPTER_GENOMIC_POSITION': accepter_genome_position
            }
            record = seqvar.VariantRecord(
                location=location,
                ref=ref_seq,
                alt='<FUSION>',
                _type='Fusion',
                _id=fusion_id,
                attrs=attrs
            )
            records.append(record)
        return records
