""" Moduel for CircRNA """
from __future__ import annotations
from typing import Iterable, List
from moPepGen.SeqFeature import SeqFeature, FeatureLocation
from moPepGen import dna, gtf


def parse(path:str) -> Iterable[CircRNAModel]:
    """ Parse a circRNA TSV file and returns an iterable of CircRNAModel
    object """
    with open(path, 'r') as handle:
        i = 0
        line = next(handle, None)

        while line:
            i += 1
            if line:
                if line.startswith('#'):
                    line = next(handle, None)
                    continue
                fields = line.rstrip().split('\t')

            gene_id = fields[0]
            start = int(fields[1])

            positions = fields[2].split(',')
            lengths = fields[3].split(',')

            if fields[4] == '.':
                introns = []
            else:
                introns = [int(x) for x in fields[4].split(',')]

            circ_id = fields[5]
            transcript_ids = fields[6].split(',')
            gene_name = fields[7]

            if len(positions) != len(lengths):
                raise ValueError(
                    f'Number of postions and lengths not match in line {i}'
                )

            fragments:List[SeqFeature] = []
            for j, (position, length) in enumerate(zip(positions, lengths)):
                position = int(position)
                length = int(length)
                start_j = start + position
                end_j = start_j + length
                location = FeatureLocation(seqname=gene_id, start=start_j, end=end_j)
                frag = SeqFeature(
                    chrom=gene_id, location=location, attributes={},
                    type='intron' if j+1 in introns else 'exon'
                )
                fragments.append(frag)

            yield CircRNAModel(gene_id, fragments, introns, circ_id,
                transcript_ids, gene_name)
            line = next(handle, None)


class CircRNAModel():
    """
    Attributes:
        gene_id (str)
        fragments (List[SeqFeature])
        intron (List[int])
        id (str)
        transcript_ids (str)
        gene_name (str)
        gene_locations (List[SeqFeature])
    """
    def __init__(self, gene_id:str, fragments:List[SeqFeature], intron:List[int],
            _id:str, transcript_ids:str,  gene_name:str):
        """ Constructor """
        self.gene_id = gene_id
        self.fragments = fragments
        self.intron = intron
        self.id = _id
        self.transcript_ids = transcript_ids
        self.gene_name = gene_name
        self.gene_locations = []

    def get_gene_coordinates(self, gene:gtf.GeneAnnotationModel) -> None:
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

    def get_circ_rna_sequence(self, seq:dna.DNASeqRecordWithCoordinates):
        """ Get the DNA sequence of the circRNA.

        Args:
            seq (dna.DNASeqRecordWithCoordinates): The DNA sequence of the
                transcript where the circRNA comes from.
        """
        circ = None
        for fragment in self.fragments:
            new_seq = seq[int(fragment.location.start):int(fragment.location.end)]
            circ = circ + new_seq if circ else new_seq
        return circ
