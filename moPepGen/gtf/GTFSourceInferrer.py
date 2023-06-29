""" Infer GTF source (e.g. GENCODE/ENSEMBL) """
from __future__ import annotations
from typing import Dict, TYPE_CHECKING


if TYPE_CHECKING:
    from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature

class GTFSourceInferrer():
    """ Infer GTF source (e.g. GENOCDE/ENSEMBL) """
    def __init__(self):
        """ Constructor """
        self.max_iter = 100
        self.data:Dict[str,int] = {}
        self.count = 0
        self.source:str = None

    def infer(self, record:GTFSeqFeature) -> str:
        """ Infer the source of a GTF record """
        if self.count > self.max_iter:
            if not self.source:
                self.source = sorted(self.data.items(), key=lambda x:x[1])[-1][0]
            return self.source
        self.count += 1
        record.infer_annotation_source()
        source = record.source
        if source not in self.data:
            self.data[source] = 1
        else:
            self.data[source] += 1
        return source
