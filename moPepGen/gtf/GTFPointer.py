""" """
from __future__ import annotations
from typing import Dict, Set, Iterable, Mapping, Union, IO, TYPE_CHECKING
from collections import deque
from moPepGen.gtf import GtfIO
from moPepGen.gtf.GeneAnnotationModel import GeneAnnotationModel
from moPepGen.gtf.TranscriptAnnotationModel import TranscriptAnnotationModel
from moPepGen.gtf.TranscriptAnnotationModel import (
    TranscriptAnnotationModel,
    GTF_FEATURE_TYPES
)

if TYPE_CHECKING:
    from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature

GENE_DICT_CACHE_SIZE = 10
TX_DICT_CACHE_SIZE = 10

class GTFPointer():
    """ """
    def __init__(self, handle:IO, key:str, start:int, end:int, source:str):
        """ """
        self.handle = handle
        self.key = key
        self.start = start
        self.end = end
        self.source = source

    def __len__(self) -> int:
        """ length """
        return self.end - self.start

class GenePointer(GTFPointer):
    """ """
    def __init__(self, handle: IO, key: str, start: int, end: int, source: str,
            transcripts:Set[str]=None):
        super().__init__(handle, key, start, end, source)
        self.transcripts = transcripts or set()

    def load(self) -> GeneAnnotationModel:
        """ """
        cur = self.handle.tell()
        offset = self.start - cur
        self.handle.seek(offset, 1)
        buffer:bytes = self.handle.read(len(self))
        lines = buffer.decode('utf-8').rstrip().split('\n')
        if len(lines) > 1:
            raise ValueError(f"Multiple lines found for gene {self.key}")
        record = GtfIO.line_to_seq_feature(lines[0])
        record.__class__ = GeneAnnotationModel
        record.exons = []
        record.transcripts = list(self.transcripts)
        return record

    def to_line(self) -> str:
        """ """
        fields = [
            self.key,
            str(self.start),
            str(self.end),
            ','.join(self.transcripts)
        ]
        return '\t'.join(fields)

class TranscriptPointer(GTFPointer):
    """ """
    def __init__(self, handle: IO, key: str, start: int, end: int, source: str,
            is_protein_coding:bool=None):
        super().__init__(handle, key, start, end, source)
        self.is_protein_coding = is_protein_coding

    def load(self) -> TranscriptAnnotationModel:
        """ """
        cur = self.handle.tell()
        offset = self.start - cur
        self.handle.seek(offset, 1)
        buffer:bytes = self.handle.read(len(self))
        lines = buffer.decode('utf-8').rstrip().split('\n')
        tx_model = TranscriptAnnotationModel()

        for line in lines:
            record = GtfIO.line_to_seq_feature(line)
            feature = record.type.lower()
            if feature not in GTF_FEATURE_TYPES:
                continue
            tx_id = record.transcript_id
            if tx_id != self.key:
                raise ValueError("Loaded transcript do not match.")
            record.id = tx_id
            tx_model.add_record(feature, record)
        tx_model.sort_records()
        return tx_model

    def to_line(self) -> str:
        """ """
        fields = [
            self.key,
            str(self.start),
            str(self.end),
            str(self.is_protein_coding)
        ]
        return "\t".join(fields)

def iterate_pointer(handle:IO, source:str=None) -> Iterable[Union[GenePointer, TranscriptPointer]]:
    """ """
    if not source:
        count = 0
        inferred = {}

    cur_gene_id:str = None
    cur_tx_id:str = None
    cur_gene_pointer = None
    cur_tx_pointer = None
    line_end = 0
    for line in handle:
        line:bytes
        line_start = line_end
        line_end += len(line)
        line = line.decode('utf-8')

        if line.startswith('#'):
            continue

        record = GtfIO.line_to_seq_feature(line)

        if not source:
            if count > 100:
                inferred = sorted(inferred.items(), key=lambda x: x[1])
                source = inferred[-1][0]
                record.source = source
            else:
                count += 1
                record.infer_annotation_source()
                inferred_source = record.source
                if inferred_source not in inferred:
                    inferred[inferred_source] = 1
                else:
                    inferred[inferred_source] += 1
        else:
            record.source = source

        if record.type.lower() == 'gene':
            if cur_gene_pointer:
                yield cur_gene_pointer
            if cur_tx_pointer:
                yield cur_tx_pointer
            cur_gene_id = record.gene_id
            cur_gene_pointer = GenePointer(handle, cur_gene_id, line_start, line_end, record.source)
            cur_tx_id = None
            cur_tx_pointer = None
        else:
            if cur_tx_id != record.transcript_id:
                if cur_tx_pointer:
                    yield cur_tx_pointer
                cur_tx_id = record.transcript_id
                cur_tx_pointer = TranscriptPointer(handle, cur_tx_id, line_start, line_end, record.source)
            else:
                cur_tx_pointer.end = line_end
            if cur_gene_pointer:
                cur_gene_pointer.transcripts.add(cur_tx_id)
    if cur_gene_pointer:
        yield cur_gene_pointer
    if cur_tx_pointer:
        yield cur_tx_pointer

class GTFPointerDict(dict, Mapping[str, GTFPointer]):
    """ """
    def __init__(self, *args:Dict[str,GTFPointer], **kwargs:Dict[str,GTFPointer]):
        """ """
        for val in kwargs.values():
            self.validate(val, GTFPointer)
        if args:
            if not isinstance(args[0], dict):
                raise TypeError('Data type invalid')
            for val in args[0].values():
                self.validate(val, GTFPointer)
        super().__init__(*args, **kwargs)

    def __iter__(self) -> Iterable[str]:
        """"""
        for k in self.keys():
            yield k

    def validate(self, val, type):
        """ """
        if not type is GTFPointer or not issubclass(type, GTFPointer):
            raise TypeError(f"type = {type} is not supported")
        if not isinstance(val, type):
            raise TypeError(
                f"{self.__class__} only accecpts {type.__class__} objects"
            )

    def __setitem__(self, __key: str, __value: GTFPointer) -> None:
        self.validate(__value, GTFPointer)
        return super().__setitem__(__key, __value)

    def get_pointer(self, __key: str) -> GTFPointer:
        """ """
        return super().__getitem__(__key)

class GenePointerDict(GTFPointerDict):
    """ """
    def __init__(self, *args: Dict[str, GenePointer], **kwargs: Dict[str, GenePointer]):
        for val in kwargs.values():
            self.validate(val, GenePointer)
        if args:
            if not isinstance(args[0], dict):
                raise TypeError('Data type invalid')
            for val in args[0].values():
                self.validate(val, GenePointer)
        super().__init__(*args, **kwargs)
        self._cache = {}
        self._cached_keys = deque()

    def __getitem__(self, __key: str) -> GeneAnnotationModel:
        if __key in self._cache:
            return self._cache[__key]
        self._cached_keys.appendleft(__key)
        if len(self._cached_keys) > GENE_DICT_CACHE_SIZE:
            key_pop = self._cached_keys.pop()
            self._cache.pop(key_pop)
        pointer:GenePointer = self.get_pointer(__key)
        val = pointer.load()
        self._cache[__key] = val
        return val

class TranscriptPointerDict(GTFPointerDict):
    """ """
    def __init__(self, *args: Dict[str, TranscriptPointer], **kwargs: Dict[str, TranscriptPointer]):
        for val in kwargs.values():
            self.validate(val, TranscriptPointer)
        if args:
            if not isinstance(args[0], dict):
                raise TypeError('Data type invalid')
            for val in args[0].values():
                self.validate(val, TranscriptPointer)
        super().__init__(*args, **kwargs)
        self._cache = {}
        self._cached_keys = deque()

    def __getitem__(self, __key: str
            ) -> TranscriptAnnotationModel:
        if __key in self._cache:
            return self._cache[__key]
        self._cached_keys.appendleft(__key)
        if len(self._cached_keys) > TX_DICT_CACHE_SIZE:
            key_pop = self._cached_keys.pop()
            self._cache.pop(key_pop)
        pointer:TranscriptPointer = self.get_pointer(__key)
        val = pointer.load()
        val.is_protein_coding = pointer.is_protein_coding is True
        val.transcript.source = pointer.source
        self._cache[__key] = val
        return val
