""" Module for GTFSeqFeature """
from __future__ import annotations
import re
from moPepGen.SeqFeature import SeqFeature


class GTFSeqFeature(SeqFeature):
    """ GTF specific SeqFeture """
    def __init__(self, *args, source:str=None, frame:int=None, **kwargs):
        """ Constructor """
        super().__init__(*args, **kwargs)
        self.source = source
        self.frame = frame

    @property
    def biotype(self) -> str:
        """ biotype """
        if self.source == 'GENCODE':
            return self.attributes['gene_type']
        if self.source == 'ENSEMBL':
            return self.attributes['gene_biotype']
        raise ValueError(f'Annotation source {self.source} not supported')

    @biotype.setter
    def biotype(self, biotype:str) -> None:
        """ biotype """
        if self.source == 'GENCODE':
            self.attributes['gene_type'] = biotype
        elif self.source == 'ENSEMBL':
            self.attributes['gene_biotype'] = biotype
        else:
            raise ValueError(f'Annotation source {self.source} not supported')

    @property
    def gene_id(self) -> str:
        """ gene ID """
        if 'gene_id' not in self.attributes:
            return None
        return self.attributes['gene_id']

    @gene_id.setter
    def gene_id(self, gene_id:str) -> None:
        """ set gene ID """
        self.attributes['gene_id'] = gene_id

    @property
    def transcript_id(self) -> str:
        """ transcript ID """
        if 'transcript_id' not in self.attributes:
            return None
        return self.attributes['transcript_id']

    @transcript_id.setter
    def transcript_id(self, transcript_id:str) -> None:
        """ set transcript ID """
        self.attributes['transcript_id'] = transcript_id

    @property
    def protein_id(self) -> str:
        """ protein ID """
        if 'protein_id' not in self.attributes:
            return None
        return self.attributes['protein_id']

    @protein_id.setter
    def protein_id(self, protein_id:str) -> None:
        """ set protein ID """
        self.attributes['protein_id'] = protein_id

    @property
    def gene_name(self) -> str:
        """ gene name """
        if 'gene_name' in self.attributes:
            return self.attributes['gene_name']
        return ''

    @gene_name.setter
    def gene_name(self, gene_name:str) -> None:
        """ set gene name """
        self.attributes['gene_name'] = gene_name

    def _shift(self, offset:int) -> GTFSeqFeature:
        """ shift by i """
        new_feature = super()._shift(offset)
        new_feature.__class__ = self.__class__
        return new_feature

    def infer_annotation_source(self):
        """ Infer the annotation source from a GTF record. Returns GENCODE,
        ENSEMBL, or None. """
        pattern = re.compile('chr[0-9XYM]{1,2}')
        if pattern.search(self.chrom):
            self.source = 'GENCODE'
            return

        pattern = re.compile(r'[A-Z]{2}[0-9]{6}\.[0-9]{1}')
        if pattern.search(self.chrom):
            self.source = 'GENCODE'
            return

        pattern = re.compile('[0-9]{1,2}')
        if pattern.search(self.chrom):
            self.source = 'ENSEMBL'
            return

        if self.chrom in ['X', 'Y', 'MT']:
            self.source = 'ENSEMBL'
            return

        pattern = re.compile('^CHR_H.+')
        if pattern.search(self.chrom):
            self.source = 'ENSEMBL'
