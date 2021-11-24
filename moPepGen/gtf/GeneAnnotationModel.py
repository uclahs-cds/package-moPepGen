""" Module for Gene Annotation Model """
from typing import List
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen import dna
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature


class GeneAnnotationModel(GTFSeqFeature):
    """ A GeneAnnotationModel object holds the genomic range of a gene, and
    the transcript IDs that it is associated with. """
    def __init__(self, chrom:str, attributes:dict, *args,
            transcripts:List[str]=None, **kwargs):
        """ Constructor """
        super().__init__(chrom=chrom, attributes=attributes, *args, **kwargs)
        if transcripts is None:
            transcripts = []
        self.transcripts = transcripts
        self.exons:List[GTFSeqFeature] = []

    def get_gene_sequence(self, chrom:dna.DNASeqRecord
            ) -> dna.DNASeqRecordWithCoordinates:
        """ Returns the DNA sequence of the gene, with all exons and introns
        reserved. """
        if self.strand == 1:
            seq = chrom.seq[self.location.start:self.location.end]
        elif self.strand == -1:
            seq = chrom.seq[self.location.start:self.location.end]
            seq = seq.reverse_complement()
        else:
            raise ValueError('Gene is unstranded.')
        gene_id = self.gene_id
        location = FeatureLocation(seqname=gene_id, start=0, end=len(seq))
        location = MatchedLocation(query=location, ref=location)
        return dna.DNASeqRecordWithCoordinates(seq=seq, locations=[location])

    def is_mrna_start_nf(self) -> bool:
        """ Whether the gene annotation has the mrna_start_NF tag """
        return 'tag' in self.attributes and \
            'mrna_start_NF' in self.attributes['tag']

    def is_cds_start_nf(self) -> bool:
        """ Whether the gene annotation has the cds_start_NF tag """
        return 'tag' in self.attributes and \
            'cds_start_NF' in self.attributes['tag']
