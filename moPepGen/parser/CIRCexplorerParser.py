""" Module for CIRCexplorer parser """
from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Tuple
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen import gtf
from moPepGen.circ import CircRNAModel


def parse(path:Path) -> Iterable[CIRCexplorerKnownRecord]:
    """ parse """
    with open(path, 'rt') as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            yield CIRCexplorerKnownRecord(
                chrom=fields[0],
                start=int(fields[1]),
                end=int(fields[2]),
                name=fields[3],
                score=float(fields[4]),
                strand=fields[5],
                thick_start=int(fields[6]),
                thick_end=int(fields[7]),
                item_rgb=[int(x) for x in fields[8].split(',')],
                exon_count=int(fields[9]),
                exon_sizes=[int(x) for x in fields[10].split(',')],
                exon_offsets=[int(x) for x in fields[11].split(',')],
                read_number=int(fields[12]),
                circ_type=fields[13],
                gene_name=fields[14],
                isoform_name=fields[15],
                index=[int(x) for x in fields[16].split(',')],
                flank_intron=fields[17]
            )


class CIRCexplorerKnownRecord():
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

    def convert_to_circ_rna(self, anno:gtf.GenomicAnnotation
            ) -> CircRNAModel:
        """ COnvert a CIRCexplorerKnownRecord to CircRNAModel. """
        tx_model = anno.transcripts[self.isoform_name]
        gene_id = tx_model.transcript.gene_id
        gene_model = anno.genes[gene_id]
        transcript_ids = [self.isoform_name]

        fragments:SeqFeature = []
        intron:List[int] = []

        if self.circ_type == 'circRNA':
            fragment_type = 'exon'
            circ_id = 'CIRC'
        elif self.circ_type == 'ciRNA':
            fragment_type = 'intron'
            circ_id = 'CI'
        else:
            raise ValueError(f'circRNA type unsupported: {self.circ_type}')

        fragment_ids = []

        for i, exon_size in enumerate(self.exon_sizes):
            exon_offset = self.exon_offsets[i]
            start = anno.coordinate_genomic_to_gene(
                self.start + exon_offset, gene_id)
            end = anno.coordinate_genomic_to_gene(
                self.start + exon_offset + exon_size - 1, gene_id)

            if gene_model.strand == -1:
                start, end = end, start
            end += 1

            location = FeatureLocation(seqname=gene_id, start=start, end=end)

            fragment = SeqFeature(
                chrom=gene_id, location=location, attributes={},
                type=fragment_type
            )

            if fragment_type == 'exon':
                exon_index = anno.find_exon_index(gene_id, fragment)
                fragment_ids.append( f"E{exon_index + 1}")

            elif fragment_type == 'intron':
                intron_index = anno.find_intron_index(gene_id, fragment)
                fragment_ids.append(f"I{intron_index + 1}")
                intron.append(i)

            fragments.append(fragment)

        fragment_ids.sort()
        circ_id += f'-{gene_id}-' + '-'.join(fragment_ids)

        return CircRNAModel(gene_id, fragments, intron, circ_id,
            transcript_ids, gene_model.gene_name)
