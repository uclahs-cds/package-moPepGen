""" Module for GTFSeqFeature """
from __future__ import annotations
import re
from moPepGen.SeqFeature import SeqFeature


class GTFSeqFeature(SeqFeature):
    """ GTF specific SeqFeture """
    def __init__(self, *args, source:str=None, **kwargs):
        """ Constructor """
        super().__init__(*args, **kwargs)
        self.source = source

    @property
    def biotype(self) -> str:
        """ biotype """
        if self.source == 'GENCODE':
            return self.attributes['gene_type']
        if self.source == 'ENSEMBL':
            return self.attributes['gene_biotype']
        raise ValueError(f'Annotation source {self.source} not supported')

    @property
    def gene_id(self) -> str:
        """ gene ID """
        return self.attributes['gene_id']

    @property
    def transcript_id(self) -> str:
        """ transcript ID """
        return self.attributes['transcript_id']

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
