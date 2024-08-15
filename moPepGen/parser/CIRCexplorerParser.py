""" Module for CIRCexplorer parser """
from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Tuple, Union
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen import gtf
from moPepGen.circ import CircRNAModel


class CIRCexplorer2KnownRecord():
    """ CIRCexplorer Known Record
    NOTE: CIRCexplorer uses 0-based and end exclusive coordinates.
    """
    def __init__(self, chrom:str, start:int, end:int, name:str, score:float,
            strand:str, thick_start:int, thick_end:int,
            item_rgb:Tuple[int,int,int], exon_count:int, exon_sizes:List[int],
            exon_offsets:List[int], read_number: int, circ_type:str,
            gene_name:str, isoform_name:str, index:List[int], flank_intron:str):
        """ Constructor """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.thick_start = thick_start
        self.thick_end = thick_end
        self.item_rgb = item_rgb
        self.exon_count = exon_count
        self.exon_sizes = exon_sizes
        self.exon_offsets = exon_offsets
        self.read_number = read_number
        self.circ_type = circ_type
        self.gene_name = gene_name
        self.isoform_name = isoform_name
        self.index = index
        self.flank_intron = flank_intron

    def convert_to_circ_rna(self, anno:gtf.GenomicAnnotation,
            intron_start_range:Tuple[int,int]=(0,0),
            intron_end_range:Tuple[int,int]=(0,0)) -> CircRNAModel:
        """ Convert a CIRCexplorerKnownRecord to CircRNAModel. """
        tx_id = self.isoform_name
        tx_model = anno.transcripts[tx_id]
        gene_id = tx_model.transcript.gene_id
        strand = tx_model.transcript.strand

        fragments:List[SeqFeature] = []
        intron:List[int] = []

        if self.circ_type == 'circRNA':
            fragment_type = 'exon'
        elif self.circ_type == 'ciRNA':
            fragment_type = 'intron'
        else:
            raise ValueError(f'circRNA type unsupported: {self.circ_type}')

        fragment_ids:List[Tuple[SeqFeature, str, int]] = []

        for i, exon_size in enumerate(self.exon_sizes):
            exon_offset = self.exon_offsets[i]
            start = anno.coordinate_genomic_to_gene(
                self.start + exon_offset, gene_id)
            end = anno.coordinate_genomic_to_gene(
                self.start + exon_offset + exon_size - 1, gene_id)

            if strand == -1:
                start, end = end, start
            end += 1

            location = FeatureLocation(
                seqname=gene_id, start=start, end=end, strand=strand
            )

            fragment = SeqFeature(
                chrom=tx_id, location=location, attributes={},
                type=fragment_type
            )

            if fragment_type == 'exon':
                exon_index = anno.find_exon_index(tx_id, fragment)
                fragment_ids.append((fragment, 'E', exon_index))

            elif fragment_type == 'intron':
                intron_index = anno.find_intron_index(
                    tx_id, fragment,
                    intron_start_range=intron_start_range,
                    intron_end_range=intron_end_range
                )
                fragment_ids.append((fragment, 'I', intron_index))
                intron.append(i)

            fragments.append(fragment)

        fragment_ids.sort(key=lambda x: x[0])

        genomic_location = f"{self.chrom}:{self.start}:{self.end}"
        start_gene = anno.coordinate_genomic_to_gene(self.start, gene_id)
        end_gene = anno.coordinate_genomic_to_gene(self.end - 1, gene_id)
        if strand == -1:
            start_gene, end_gene = end_gene, start_gene
        end_gene += 1
        backsplicing_site = FeatureLocation(
            seqname=tx_model.transcript.chrom,
            start=start_gene,
            end=end_gene
        )
        circ_id = f"CIRC-{tx_id}-{start_gene}:{end_gene}"

        return CircRNAModel(
            transcript_id=tx_id,
            fragments=fragments,
            intron=intron,
            _id=circ_id,
            gene_id=gene_id,
            gene_name=tx_model.transcript.gene_name,
            genomic_location=genomic_location,
            backsplicing_site=backsplicing_site
        )

    def is_valid(self, min_read_number:int) -> bool:
        """ Check if this is valid CIRCexplorer record """
        return self.read_number >= min_read_number

class CIRCexplorer3KnownRecord(CIRCexplorer2KnownRecord):
    """ CIRCexplorer3 """
    def __init__(self, *args, fpb_circ:float, fpb_linear:float,
            circ_score:float, **kwargs):
        super().__init__(*args, **kwargs)
        self.fpb_circ = fpb_circ
        self.fpb_linear = fpb_linear
        self.circ_score = circ_score

    def is_valid(self, min_read_number: int, min_fbr_circ:float=None,
            min_circ_score:float=None) -> bool:
        """ Check if is valid """
        #pylint: disable=W0221
        if min_fbr_circ and self.fpb_circ < min_fbr_circ:
            return False
        if min_circ_score and self.circ_score < min_circ_score:
            return False
        return super().is_valid(min_read_number)


CIRCexplorerKnownRecord = Union[CIRCexplorer2KnownRecord, CIRCexplorer3KnownRecord]

def parse(path:Path, circ_explorer_3:bool=False
        ) -> Iterable[CIRCexplorerKnownRecord]:
    """ parse """
    with open(path, 'rt') as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            data = {
                'chrom': fields[0],
                'start': int(fields[1]),
                'end': int(fields[2]),
                'name': fields[3],
                'score': float(fields[4]),
                'strand': fields[5],
                'thick_start': int(fields[6]),
                'thick_end': int(fields[7]),
                'item_rgb': [int(x) for x in fields[8].split(',')],
                'exon_count': int(fields[9]),
                'exon_sizes': [int(x) for x in fields[10].split(',')],
                'exon_offsets': [int(x) for x in fields[11].split(',')],
                'read_number': int(fields[12]),
                'circ_type': fields[13],
                'gene_name': fields[14],
                'isoform_name': fields[15],
                'index': [int(x) for x in fields[16].split(',')],
                'flank_intron': fields[17]
            }
            if circ_explorer_3:
                yield CIRCexplorer3KnownRecord(
                    fpb_circ=float(fields[18]),
                    fpb_linear=float(fields[19]),
                    circ_score=float(fields[20]),
                    **data
                )
            else:
                yield CIRCexplorer2KnownRecord(**data)
