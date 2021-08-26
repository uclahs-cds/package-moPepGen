""" Moduel for CircRNA """
from __future__ import annotations
from typing import List, TYPE_CHECKING
from moPepGen.SeqFeature import SeqFeature, FeatureLocation


if TYPE_CHECKING:
    from moPepGen.dna import DNASeqRecordWithCoordinates
    from moPepGen.gtf import GeneAnnotationModel

class CircRNAModel():
    """
    Attributes:
        gene_id (str)
        fragments (List[SeqFeature])
        intron (List[int])
        id (str)
        transcript_ids (List[str])
        gene_name (str)
        gene_locations (List[SeqFeature])
    """
    def __init__(self, gene_id:str, fragments:List[SeqFeature], intron:List[int],
            _id:str, transcript_ids:List[str],  gene_name:str):
        """ Constructor """
        self.gene_id = gene_id
        self.fragments = fragments
        self.intron = intron
        self.id = _id
        self.transcript_ids = transcript_ids
        self.gene_name = gene_name
        self.gene_locations = []

    def get_gene_coordinates(self, gene:GeneAnnotationModel) -> None:
        """ Get the coordinates of the gene """
        features:List[SeqFeature] = []
        for i,fragment in enumerate(self.fragments):
            start = fragment.location.start
            end = fragment.location.end
            feature = SeqFeature(
                chrom=gene.id,
                location=FeatureLocation(seqname=gene.id, start=start, end=end),
                attributes={},
                type='intron' if i + 1 in self.intron else 'exon',
                strand=gene.strand
            )
            features.append(feature)
        self.gene_locations = features

    def get_circ_rna_sequence(self, seq:DNASeqRecordWithCoordinates):
        """ Get the DNA sequence of the circRNA.

        Args:
            seq (DNASeqRecordWithCoordinates): The DNA sequence of the
                transcript where the circRNA comes from.
        """
        circ = None
        for fragment in self.fragments:
            new_seq = seq[int(fragment.location.start):int(fragment.location.end)]
            circ = circ + new_seq if circ else new_seq
        return circ

    def to_string(self) -> str:
        """ Convert to a string """
        gene_id = self.gene_id
        start = int(self.fragments[0].location.start)
        offset, length = [], []
        for fragment in self.fragments:
            offset.append(str(fragment.location.start - start))
            length.append(str(fragment.location.end - fragment.location.start))
        start = str(start)
        offset = ','.join(offset)
        length = ','.join(length)
        intron = ','.join([str(x) for x in self.intron])
        circ_id = self.id
        tx_id = ','.join(self.transcript_ids)
        gene_name = self.gene_name
        return f'{gene_id}\t{start}\t{offset}\t{length}\t{intron}\t{circ_id}'\
            + f'\t{tx_id}\t{gene_name}'
