""""""
from __future__ import annotations
from moPepGen.seqvar.VariantRecord import VariantRecord
from typing import List, Tuple, Iterable
from moPepGen import seqvar, dna


class REDIToolsRecord():
    """"""
    def __init__(self, region:str, position:int, reference:str, strand:int,
            coverage_q30:int, mean_quality:float, base_count:List[int],
            all_subs:List[Tuple(str, str)], frequency:float):
        self.region = region
        self.position = position
        self.reference = reference
        self.starnd = strand
        self.all_subs = all_subs
        self.coverage_q30 = coverage_q30
        self.mean_quality = mean_quality
        self.base_count = base_count
        self.all_subs = all_subs
        self.frequency = frequency

    @staticmethod
    def parse(path:str, mode:str='r') -> Iterable[REDIToolsRecord]:
        """"""
        with open(path, mode) as handle:
            next(handle)
            for line in handle:
                fields = line.rstrip().split('\t')
                base_count = [int(i) for i in fields[6].strip('][')\
                        .split(', ')]
                all_subs = []
                for sub in fields[7].split(' '):
                    if len(sub) > 2:
                        raise ValueError('length of sub larger than 2')
                    all_subs.append((sub[0], sub[1]))
                yield REDIToolsRecord(
                    region=fields[0],
                    position=int(fields[1]),
                    reference=fields[2],
                    strand=int(fields[3]),
                    coverage_q30=int(fields[4]),
                    mean_quality=float(fields[5]),
                    base_count=base_count,
                    all_subs=all_subs,
                    frequency=float(fields[8])
                )

    def convert_to_variant_record(self, seq:dna.DNASeqRecord
            ) -> seqvar.VariantRecord:
        """"""
        # need annotation
        # potential solution: pybedtools
        VariantRecord()