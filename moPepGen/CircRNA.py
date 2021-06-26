""" Moduel for CircRNA """
from __future__ import annotations
from typing import Iterable, List
from moPepGen.SeqFeature import SeqFeature, FeatureLocation
from moPepGen import dna


def parse(path:str) -> Iterable[CircRNAModel]:
    """ Parse a circRNA BED file and returns an iterable of CircRNAModel
    object """
    with open(path, 'r') as handle:
        line = next(handle, None)
        transcript_id = None
        circ_rna_id = None
        circ = None

        while True:
            if line:
                if line.startswith('#'):
                    line = next(handle, None)
                    continue
                fields = line.rstrip().split('\t')

            if transcript_id and ( not line or (transcript_id != fields[0]\
                    and circ_rna_id != fields[3])):
                yield circ
                if not line:
                    return
                transcript_id = None

            location = FeatureLocation(seqname=fields[0], start=int(fields[1]),
                end=int(fields[2]))
            exon = SeqFeature(chrom=fields[0], location=location, attributes={})

            if transcript_id is None:
                transcript_id = fields[0]
                circ_rna_id = fields[3]
                circ = CircRNAModel(transcript_id=transcript_id, exons=[exon],
                    _id=fields[3], gene_id=fields[4], gene_name=fields[5])
            elif transcript_id == fields[0]:
                circ.exons.append(exon)
            line = next(handle, None)


class CircRNAModel():
    """
    Attributes:
        transcript_id (str)
        exons (List[SeqFeature])
        id (str)
        gene_id (str)
        gene_name (str)
    """
    def __init__(self, transcript_id:str, exons:List[SeqFeature],
            _id:str, gene_id:str, gene_name:str):
        """ Constructor """
        self.transcript_id = transcript_id
        self.exons = exons
        self.id = _id
        self.gene_id = gene_id
        self.gene_name = gene_name

    def get_circ_rna_sequence(self, seq:dna.DNASeqRecordWithCoordinates):
        """ Get the DNA sequence of the circRNA.

        Args:
            seq (dna.DNASeqRecordWithCoordinates): The DNA sequence of the
                transcript where the circRNA comes from.
        """
        circ = None
        for exon in self.exons:
            new_seq = seq[int(exon.location.start):int(exon.location.end)]
            circ = circ + new_seq if circ else new_seq
        return circ
