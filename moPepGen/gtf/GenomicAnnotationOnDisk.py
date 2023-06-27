""" Module for GenomicAnnotationOnDisk """
from __future__ import annotations
from typing import Dict, IO, Union, Tuple
import io
import gzip
import lzma
from pathlib import Path
from moPepGen.version import MetaVersion
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation
from moPepGen.gtf.GTFPointer import (
    GenePointerDict, TranscriptPointerDict, iterate_pointer,
    GenePointer, TranscriptPointer
)


class GenomicAnnotationOnDisk(GenomicAnnotation):
    """ """
    def __init__(self, genes:GenePointerDict=None,
            transcripts:TranscriptPointerDict=None,
            source:str=None):
        """ """
        self.genes = genes or GenePointerDict()
        self.transcripts = transcripts or TranscriptPointerDict()
        self.source = source
        self.gene_id_version_mapper:Dict[str,str] = None
        self.version = MetaVersion()
        self._cached_tx_seqs = []
        self.handle = None

    def __del__(self):
        """ """
        if self.handle:
            self.handle.close()

    def dump_gtf(self, *args, **kwargs) -> None:
        raise NotImplementedError

    def infer_source(self):
        """ """
        i = 0
        inferred = {}
        for tx in self.transcripts.keys():
            i += 1
            if i > 100:
                break
            s = self.transcripts.get_pointer(tx).source
            if s in inferred:
                inferred[s] += 1
            else:
                inferred[s] = 1
        self.source = sorted(inferred.items(), key=lambda it: it[1])[-1][0]

    def init_handle(self, handle:Union[str, IO, Path]):
        """ """
        if isinstance(handle, str):
            handle = Path(handle)

        if isinstance(handle, Path):
            if handle.suffix.lower() == '.gtf':
                ihandle = open(handle, 'rb')
            elif handle.suffix.lower() == '.lz':
                ihandle = lzma.open(handle, 'rb')
            elif handle.suffix.lower() == '.gz':
                ihandle = gzip.open(handle, 'rb')
        elif isinstance(handle, io.IOBase):
            ihandle = handle
        else:
            raise ValueError(f"handle {handle} of type {handle.__class__} is not supported.")

        self.handle = ihandle

    @staticmethod
    def get_index_files(file:Union[str,Path]) -> Tuple[Path,Path]:
        """ """
        if isinstance(file, str):
            file = Path(file)
        basename = str(file.name).split('.gtf')[0]
        gene_idx_file = file.parent/f"{basename}_gene.idx"
        tx_idx_file = file.parent/f"{basename}_tx.idx"
        return gene_idx_file, tx_idx_file

    def generate_index(self, handle:Union[IO, str, Path], source:str=None):
        """ """
        self.init_handle(handle)

        for pointer in iterate_pointer(self.handle, source):
            if isinstance(pointer, GenePointer):
                self.genes[pointer.key] = pointer
            else:
                self.transcripts[pointer.key] = pointer

        if source:
            self.source = source
        else:
            self.infer_source()

    def load_index(self, file:Union[str,Path]):
        """ """
        self.init_handle(file)
        gene_idx_file, tx_idx_file = self.get_index_files(file)

        gene_source = None
        with open(gene_idx_file, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    line = line.lstrip('#').lstrip()
                    if line.startswith('source='):
                        gene_source = line.split('=')[1]
                    continue

                if not gene_source:
                    raise ValueError("Cannot find source from gene idx.")

                fields = line.rstrip().split('\t')
                pointer = GenePointer(
                    self.handle, key=fields[0],
                    start=int(fields[1]), end=int(fields[2]),
                    source=gene_source, transcripts=fields[3].split(',')
                )
                self.genes[pointer.key] = pointer

        tx_source = None
        with open(tx_idx_file, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    line = line.lstrip('#').lstrip()
                    if line.startswith('source='):
                        tx_source = line.split('=')[1]
                    continue

                if not tx_source:
                    raise ValueError("Cannot find source from gene idx.")

                if tx_source != gene_source:
                    raise ValueError('tx_source does not match with gene_source')

                fields = line.rstrip().split('\t')
                is_protein_coding = {'None': None, 'True': True, 'False': False}[fields[3]]
                pointer = TranscriptPointer(
                    handle=self.handle, key=fields[0],
                    start=int(fields[1]), end=int(fields[2]),
                    source=tx_source, is_protein_coding=is_protein_coding
                )
                self.transcripts[pointer.key] = pointer
