""" This module defines the class logic for the GTF annotations.
"""
from typing import List, Dict
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen import seqvar
from . import GtfIO
from .TranscriptAnnotationModel import TranscriptAnnotationModel, _GTF_FEATURE_TYPES
from .GeneAnnotationModel import GeneAnnotationModel

class GenomicAnnotation():
    """ A VEPTranscripts object is a dict-like object that holds the GTF
    results, that the keys are the transcript IDs and values are instances of
    VEPRecord.
    """
    def __init__(self, genes:Dict[str,GeneAnnotationModel]=None,
            transcripts:Dict[str, TranscriptAnnotationModel]=None):
        """ Construct a GenomicAnnotation object """
        if genes is None:
            genes = {}
        if transcripts is None:
            transcripts = {}
        self.genes = genes
        self.transcripts = transcripts

    def __repr__(self) -> str:
        """ Return a string representation """
        result = ""
        i = 0
        it = iter(self.genes.keys())
        gene_id = next(it, None)
        while i < len(self.genes) and gene_id:
            result += f"'{gene_id}': {self.genes[gene_id].__repr__()}"
            if i == 3 and len(self.genes) > 7:
                result += "...\n"
                i = len(self.transcripts) - 4
            i += 1
            gene_id = next(it, None)
        result += f'\n{len(self.genes)} genes'

        i = 0
        it = iter(self.transcripts.keys())
        key = next(it, None)
        while i < len(self) and key:
            result += f"'{key}': {self.transcripts[key].__repr__()}\n"
            if i == 3 and len(self.transcripts) > 7:
                result += "...\n"
                i = len(self.transcripts) - 4
            i += 1
            key = next(it, None)
        result += f'\n{len(self.transcripts)} transcripts'
        return result

    def dump_gtf(self, path:str, biotype:List[str]=None)->None:
        """ Dump a GTF file into a GenomicAnnotation

        Args:
            path (str): Path to a GTF file.
            biotype (List[str]): The annotation biotype to keep. Features
                in a GTF can be annotated as protein_coding, miRNA, lncRNA,
                miRNA, etc.
        """
        record:SeqFeature
        for record in GtfIO.parse(path):
            if biotype is not None and record.biotype not in biotype:
                continue
            feature = record.type.lower()
            if feature == 'gene':
                gene_id = record.attributes['gene_id']
                if gene_id in self.genes:
                    raise ValueError(f'Same gene has multiple records: {gene_id}')
                record.__class__ = GeneAnnotationModel
                record.transcripts = []
                self.genes[gene_id] = record
                continue

            if feature not in _GTF_FEATURE_TYPES:
                continue
            transcript_id = record.attributes['transcript_id']
            record.id = transcript_id
            if transcript_id not in self.transcripts.keys():
                self.transcripts[transcript_id] = TranscriptAnnotationModel()
            self.transcripts[transcript_id].add_record(feature, record)

            gene_id = record.attributes['gene_id']
            if gene_id not in self.genes:
                raise ValueError(f'Gene ID {gene_id} not found')
            if transcript_id not in self.genes[gene_id].transcripts:
                self.genes[gene_id].transcripts.append(transcript_id)

        for transcript_model in self.transcripts.values():
            transcript_model.sort_records()

    def variant_coordinates_to_gene(self, variant:seqvar.VariantRecord,
            gene_id:str) -> seqvar.VariantRecord:
        """ Convert the coordinates of variant from transcript to gene """
        transcript_id = variant.location.seqname
        start = variant.location.start
        end = variant.location.end

        if transcript_id not in self.genes[gene_id].transcripts:
            raise ValueError("The variant isn't associated with the gene.")

        transcript_model = self.transcripts[transcript_id]

        if transcript_model.transcript.location.strand == 1:
            it = iter(transcript_model.exon)
            exon = next(it, None)
            left = 0
            while exon:
                right = left + exon.location.end - exon.location.start
                if right <= start:
                    left = right
                    exon = next(it, None)
                    continue
                start_genomic = start - left + exon.location.start
                break

            while exon:
                right = left + exon.location.end - exon.location.start
                if right < end:
                    left = right
                    exon = next(it, None)
                    continue
                end_genomic = end - left + exon.location.start
                break

            start_gene = start_genomic - self.genes[gene_id].location.start
            end_gene = end_genomic - self.genes[gene_id].location.start

        elif transcript_model.transcript.location.strand == -1:
            it = reversed(transcript_model)
            exon = next(it, None)
            right = 0
            while exon:
                left = right + exon.location.end - exon.location.start
                if left <= start:
                    right = left
                    exon = next(it, None)
                    continue
                start_genomic = left - start + exon.location.start
                break

            while exon:
                left = right + exon.location.end - exon.location.start
                if left <= end:
                    right = left
                    exon = next(it, None)
                    continue
                end_genomic = left - end + exon.location.start
                break

            start_gene = self.genes[gene_id].location.start - start_genomic
            end_gene = self.genes[gene_id].location.start - end_genomic

        else:
            raise ValueError('Transcript should not be unstranded.')

        if end_gene - start_gene != end - start:
            raise ValueError('Variant seems to be over a splice site.')

        return seqvar.VariantRecord(
            location=FeatureLocation(
                seqname=gene_id,
                start=start_gene,
                end=end_gene
            ),
            ref=variant.ref,
            alt=variant.alt,
            _type=variant.type,
            _id=variant.id
        )
