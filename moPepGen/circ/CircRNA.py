""" Module for CircRNA """
from __future__ import annotations
from typing import List, TYPE_CHECKING
from moPepGen.SeqFeature import SeqFeature, FeatureLocation


if TYPE_CHECKING:
    from moPepGen.dna import DNASeqRecordWithCoordinates
    from moPepGen.gtf import GeneAnnotationModel

class CircRNAModel():
    """
    Attributes:
        transcript_id (str)
        fragments (List[SeqFeature])
        intron (List[int])
        id (str)
        gene_id (str)
        gene_name (str)
        gene_locations (List[SeqFeature])
    """
    def __init__(self, transcript_id:str, fragments:List[SeqFeature],
            intron:List[int], _id:str, gene_id:str, gene_name:str,
            genomic_location:str=''):
        """ Constructor """
        self.gene_id = gene_id
        self.fragments = fragments
        self.intron = intron
        self.id = _id
        self.transcript_id = transcript_id
        self.gene_name = gene_name
        self.gene_locations = []
        self.genomic_position = genomic_location

    def get_gene_coordinates(self, gene:GeneAnnotationModel) -> None:
        """ Get the coordinates of the gene """
        features:List[SeqFeature] = []
        for i,fragment in enumerate(self.fragments):
            start = fragment.location.start
            end = fragment.location.end
            feature = SeqFeature(
                chrom=gene.id,
                location=FeatureLocation(
                    seqname=gene.id, start=start,
                    end=end, strand=gene.strand
                ),
                attributes={},
                type='intron' if i + 1 in self.intron else 'exon'
            )
            features.append(feature)
        self.gene_locations = features

    def get_circ_rna_sequence(self, seq:DNASeqRecordWithCoordinates):
        """ Get the DNA sequence of the circRNA.

        Args:
            seq (DNASeqRecordWithCoordinates): The DNA sequence of the
                gene where the circRNA comes from.
        """
        circ = None
        for fragment in sorted(self.fragments):
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
        tx_id = self.transcript_id
        gene_name = self.gene_name
        info = f'OFFSET={offset};LENGTH={length};INTRON={intron};' +\
            f'TRANSCRIPT_ID={tx_id};GENE_SYMBOL={gene_name};' +\
            f'GENOMIC_POSITION={self.genomic_position}'
        return '\t'.join([gene_id, start, circ_id, '.', '.', '.', '.', info])
