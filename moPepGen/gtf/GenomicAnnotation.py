""" This module defines the class logic for the GTF annotations.
"""
import re
from typing import List, Dict
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen import seqvar
from . import GtfIO
from .TranscriptAnnotationModel import TranscriptAnnotationModel, GTF_FEATURE_TYPES
from .GeneAnnotationModel import GeneAnnotationModel


class GenomicAnnotation():
    """ This defines the annotation of genes and transcripts of the genome,
    reading from a GTF file.

    Attributes:
        genes (Dict[str,GeneAnnotationModel]): Keys are gene IDs and values are
            gene annotation models.
        transcripts (Dict[str, TranscriptAnnotationModel]): Keys are transcript
            IDs and values are transcript annotation models.
    """
    FAILED_TO_FIND_EXON_ERROR = 'Failed to find the exon.'
    FAILED_TO_FIND_INTRON_ERROR = 'Failed to find the intron.'
    def __init__(self, genes:Dict[str,GeneAnnotationModel]=None,
            transcripts:Dict[str, TranscriptAnnotationModel]=None,
            source:str=None):
        """ Construct a GenomicAnnotation object """
        if genes is None:
            genes = {}
        if transcripts is None:
            transcripts = {}
        self.genes = genes
        self.transcripts = transcripts
        self.source = source
        self.gene_id_version_mapper:Dict[str, str] = None

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

    def add_gene_record(self, record:SeqFeature) -> None:
        """ Add a gene record """
        gene_id = record.attributes['gene_id']
        if gene_id in self.genes:
            raise ValueError(f'Same gene has multiple records: {gene_id}')
        record.__class__ = GeneAnnotationModel
        record.exons = []
        record.transcripts = []
        self.genes[gene_id] = record

    def add_transcript_record(self, record:SeqFeature) -> None:
        """ Add a transcript record """
        feature = record.type.lower()
        if feature not in GTF_FEATURE_TYPES:
            return
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

    @staticmethod
    def infer_annotation_source(record:SeqFeature) -> str:
        """ Infer the annotation source from a GTF record. Returns GENCODE,
        ENSEMBL, or None. """
        pattern = re.compile('chr[0-9XYM]{1,2}')
        if pattern.search(record.chrom):
            return 'GENCODE'

        pattern = re.compile(r'[A-Z]{2}[0-9]{6}\.[0-9]{1}')
        if pattern.search(record.chrom):
            return 'GENCODE'

        pattern = re.compile('[0-9]{1,2}')
        if pattern.search(record.chrom):
            return 'ENSEMBL'

        if record.chrom in ['X', 'Y', 'MT']:
            return 'ENSEMBL'

        pattern = re.compile('^CHR_H.+')
        if pattern.search(record.chrom):
            return 'ENSEMBL'

        return None


    def dump_gtf(self, path:str, biotype:List[str]=None, source:str=None)->None:
        """ Dump a GTF file into a GenomicAnnotation

        Args:
            path (str): Path to a GTF file.
            biotype (List[str]): The annotation biotype to keep. Features
                in a GTF can be annotated as protein_coding, miRNA, lncRNA,
                miRNA, etc.
        """
        record:SeqFeature
        if not source:
            count = 0
            inferred = {}
        for record in GtfIO.parse(path):
            if biotype is not None and record.biotype not in biotype:
                continue

            if not source and count > 100:
                inferred = sorted(inferred.items(), key=lambda x: x[1])
                source = inferred[-1][0]

            if not source:
                count += 1
                inferred_source = self.infer_annotation_source(record)
                if inferred_source not in inferred:
                    inferred[inferred_source] = 1
                else:
                    inferred[inferred_source] += 1

            feature = record.type.lower()
            if feature == 'gene':
                self.add_gene_record(record)
                continue

            self.add_transcript_record(record)

        self.source = source

        for transcript_model in self.transcripts.values():
            transcript_model.sort_records()

    def coordinate_gene_to_transcript(self, index:int, gene:str,
            transcript:str) -> int:
        """ For a given gene position, find it the corresponding transcript
        position.

        Args:
            index (int): The gene position index.
            gene (str): The gene ID.
            transcript (str): The transcript ID to map to.
        """
        genomic_index = self.coordinate_gene_to_genomic(index, gene)
        gene_model = self.genes[gene]
        if transcript not in gene_model.transcripts:
            raise ValueError(f'The transcript {transcript} is not associated'
                f'with the gene {gene}')
        transcript_model = self.transcripts[transcript]

        return transcript_model.get_transcript_index(genomic_index)

    def coordinate_genomic_to_gene(self, index:int, gene:str) -> int:
        """ Get the gene coordinate from genomic coordinate. """
        gene_location = self.genes[gene].location
        if gene_location.strand == 1:
            return index - gene_location.start
        if gene_location.strand == -1:
            return gene_location.end - 1 - index
        raise ValueError("Don't know how to handle unstranded gene.")

    def coordinate_gene_to_genomic(self, index:int, gene:str) -> int:
        """ Get the genomic coordinate from gene coordinate """
        location = self.genes[gene].location
        if location.strand == 1:
            return location.start + index
        if location.strand == -1:
            return location.end - 1 - index
        raise ValueError("Don't know how to handle unstranded gene.")

    def coordinate_transcript_to_genomic(self, index:int, transcript:str) -> int:
        """ Get the genomic coordinate from transcript coordinate """
        tx_model = self.transcripts[transcript]
        tx_size = tx_model.transcript.location.end - tx_model.transcript.location.start
        if tx_size < index:
            raise ValueError('Index out of range.')

        if tx_model.transcript.strand == 1:
            for exon in tx_model.exon:
                size = exon.location.end - exon.location.start
                if index < size:
                    return index + exon.location.start
                index -= size
        if tx_model.transcript.strand == -1:
            for exon in reversed(tx_model.exon):
                size = exon.location.end - exon.location.start
                if index < size:
                    return exon.location.end - 1 - index
                index -= size
        raise ValueError("Don't know how to handle unstranded transcript.")

    def variant_coordinates_to_gene(self, variant:seqvar.VariantRecord,
            gene_id:str) -> seqvar.VariantRecord:
        """ Convert the coordinates of variant from transcript to gene """
        transcript_id = variant.location.seqname
        start = variant.location.start
        end = variant.location.end

        if transcript_id not in self.genes[gene_id].transcripts:
            raise ValueError("The variant isn't associated with the gene.")

        transcript_model = self.transcripts[transcript_id]
        start_genomic, end_genomic = None, None

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

            if start_genomic is None or end_genomic is None:
                raise ValueError('The variant is not is the range of the gene.')

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

        if start_genomic is None or end_genomic is None:
            raise ValueError('The variant is out off the range of the '
            'transcript.')

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

    def feature_coordinate_gene_to_genomic(self, feature:SeqFeature,
            gene_id:str) -> SeqFeature:
        """ For a given feature, converts the coordinate from gene to
        genomic. """
        gene_model = self.genes[gene_id]
        start = self.coordinate_gene_to_genomic(feature.location.start, gene_id)
        end = self.coordinate_gene_to_genomic(feature.location.end, gene_id)
        strand = gene_model.location.strand
        if strand == -1:
            start, end = end, start
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=strand)
        return SeqFeature(chrom=gene_id, location=location, attributes={})

    def feature_coordinate_genomic_to_gene(self, feature:SeqFeature,
            gene_id:str) -> SeqFeature:
        """ For a given feature, converts the coordinate from genomic to gene.
        """
        gene_model = self.genes[gene_id]
        start = self.coordinate_genomic_to_gene(feature.location.start, gene_id)
        end = self.coordinate_genomic_to_gene(feature.location.end, gene_id)
        strand = gene_model.location.strand
        if strand == -1:
            start, end = end, start
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=0)
        return SeqFeature(chrom=gene_id, location=location, attributes={})

    def create_gene_id_version_mapper(self) -> None:
        """ Create the a dict that keys are the unversioned gene ID, and values
        are the versioned. """
        self.gene_id_version_mapper = {}
        for versioned in self.genes.keys():
            unversioned = versioned.split('.')[0]
            if unversioned in self.gene_id_version_mapper:
                raise ValueError('Unversioned gene ID collapsed.')
            self.gene_id_version_mapper[unversioned] = versioned

    def get_gene_model_from_unversioned_id(self, gene_id:str) -> GeneAnnotationModel:
        """ Get the gene annotation model from an unversioned gene ID. In
        general the source and version of the annotation GTF file used should
        be consistent with the one used to call genomic variants (e.g. VEP &
        fusion). But this is useful for FusionCatcher because it uses ENSEMBL
        as default.

        Args:
            gene_id (str): The unversioned gene ID (ENSEMBL's style).

        Returns:
            The GeneAnnotationModel object of the gene.
        """
        if self.source == 'ENSEMBL':
            return self.genes[gene_id]

        if self.gene_id_version_mapper is None:
            self.create_gene_id_version_mapper()

        versioned_gene_id = self.gene_id_version_mapper[gene_id]

        return self.genes[versioned_gene_id]

    def get_all_exons_of_gene(self, gene_id) -> List[SeqFeature]:
        """ Get all exons of a gene """
        exons = set()
        gene_model = self.genes[gene_id]
        for tx_id in gene_model.transcripts:
            tx_model = self.transcripts[tx_id]
            exons.update(tx_model.exon)
        exons = list(exons)
        exons.sort()
        gene_model.exons = exons

    def find_exon_index(self, gene_id:str, feature:SeqFeature,
            coordinate:str='gene') -> int:
        """ Find the exon index of the gene.

        Args:
            gene_id (str): The gene ID.
            feature (SeqFeature): The exon to look up.
            coordinate (str): The coordinate type of the given feature. Can be
                one of 'genomic' or 'gene'. Defaults to gene.
        """
        gene_model = self.genes[gene_id]
        if not gene_model.exons:
            self.get_all_exons_of_gene(gene_id)
        exons = gene_model.exons

        strand = self.genes[gene_id].strand

        if coordinate == 'gene':
            feature = self.feature_coordinate_gene_to_genomic(feature, gene_id)
        elif coordinate != 'genomic':
            raise ValueError("Don't know how to handle the coordinate.")

        if strand == 1:
            for i, exon in enumerate(exons):
                if exon == feature:
                    return i
                if exon > feature:
                    break
        elif strand == -1:
            for i, exon in enumerate(reversed(exons)):
                if exon == feature:
                    return i
                if exon < feature:
                    break
        raise ValueError(self.FAILED_TO_FIND_EXON_ERROR)

    def find_intron_index(self, gene_id:str, feature:SeqFeature,
            coordinate:str="gene") -> int:
        """ Find the intron index of the gene.

        Args:
            gene_id (str): The gene ID.
            feature (SeqFeature): The intron to look up.
            coordinate (str): The coordinate type of the given feature. Can be
                one of 'genomic' or 'gene'. Defaults to gene.
        """
        gene_model = self.genes[gene_id]
        if not gene_model.exons:
            self.get_all_exons_of_gene(gene_id)
        exons = gene_model.exons

        strand = self.genes[gene_id].strand

        if coordinate == 'gene':
            feature = self.feature_coordinate_gene_to_genomic(feature, gene_id)

        elif coordinate != 'genomic':
            raise ValueError("Don't know how to handle the coordinate.")

        if strand == 1:
            i = 0
            while i < len(exons):
                exon = exons[i]
                if exon.location.end == feature.location.start:
                    i += 1
                    if i >= len(exons):
                        break
                    exon = exons[i]
                    if exon.location.start == feature.location.end:
                        return i - 1
                    break
                if exon.location.start > feature.location.end:
                    break
                i += 1
        elif strand == -1:
            i = 0
            j = len(exons) - i - 1
            while i < len(exons):
                exon = exons[j]
                if exon.location.start == feature.location.end:
                    i += 1
                    if i >= len(exons):
                        break
                    j = len(exons) - i - 1
                    exon = exons[j]
                    if exon.location.end == feature.location.start:
                        return i - 1
                    break
                if exon.location.end < feature.location.start:
                    break
                i += 1
                j = len(exons) - i - 1
        raise ValueError(self.FAILED_TO_FIND_INTRON_ERROR)
