""" This module defines the class logic for the GTF annotations.
"""
from __future__ import annotations
from typing import List, Dict, Tuple, TYPE_CHECKING
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen import seqvar, err
from moPepGen.version import MetaVersion
from . import GtfIO
from .TranscriptAnnotationModel import TranscriptAnnotationModel, GTF_FEATURE_TYPES
from .GeneAnnotationModel import GeneAnnotationModel
from .GTFSeqFeature import GTFSeqFeature


if TYPE_CHECKING:
    from moPepGen.aa import AminoAcidSeqDict

class GenomicAnnotation():
    """ This defines the annotation of genes and transcripts of the genome,
    reading from a GTF file.

    Attributes:
        genes (Dict[str,GeneAnnotationModel]): Keys are gene IDs and values are
            gene annotation models.
        transcripts (Dict[str, TranscriptAnnotationModel]): Keys are transcript
            IDs and values are transcript annotation models.
        source (str): Source of the genomic annotation. E.g., ENSEMBL or
            GENCODE.
    """
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
        self.version = MetaVersion()

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

    def add_gene_record(self, record:GTFSeqFeature) -> None:
        """ Add a gene record """
        gene_id = record.gene_id
        if gene_id in self.genes:
            raise ValueError(f'Same gene has multiple records: {gene_id}')
        record.__class__ = GeneAnnotationModel
        record.exons = []
        record.transcripts = []
        self.genes[gene_id] = record

    def add_transcript_record(self, record:GTFSeqFeature) -> None:
        """ Add a transcript record """
        feature = record.type.lower()
        if feature not in GTF_FEATURE_TYPES:
            return
        transcript_id = record.transcript_id
        record.id = transcript_id
        if transcript_id not in self.transcripts.keys():
            self.transcripts[transcript_id] = TranscriptAnnotationModel()
        self.transcripts[transcript_id].add_record(feature, record)

        gene_id = record.gene_id
        if gene_id not in self.genes:
            raise ValueError(f'Gene ID {gene_id} not found')
        if transcript_id not in self.genes[gene_id].transcripts:
            self.genes[gene_id].transcripts.append(transcript_id)

    def dump_gtf(self, path:str, biotype:List[str]=None, source:str=None)->None:
        """ Dump a GTF file into a GenomicAnnotation

        Args:
            path (str): Path to a GTF file.
            biotype (List[str]): The annotation biotype to keep. Features
                in a GTF can be annotated as protein_coding, miRNA, lncRNA,
                miRNA, etc.
        """
        record:GTFSeqFeature
        if not source:
            count = 0
            inferred = {}
        for record in GtfIO.parse(path):
            if biotype is not None and record.biotype not in biotype:
                continue

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

            feature = record.type.lower()
            if feature == 'gene':
                self.add_gene_record(record)
                continue

            self.add_transcript_record(record)

        if not source:
            inferred = sorted(inferred.items(), key=lambda x: x[1])
            source = inferred[-1][0]

        self.source = source

        for transcript_model in self.transcripts.values():
            transcript_model.sort_records()

    def check_protein_coding(self, proteome:AminoAcidSeqDict) -> None:
        """ Checks if each transcript is protein coding """
        for tx_id, tx_model in self.transcripts.items():
            tx_model.is_protein_coding = tx_id in proteome

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
        start_genomic = self.coordinate_transcript_to_genomic(start, transcript_id)
        end_genomic = self.coordinate_transcript_to_genomic(end - 1, transcript_id)
        if transcript_model.transcript.strand == -1:
            start_genomic, end_genomic = end_genomic, start_genomic
        end_genomic += 1
        start_gene = self.coordinate_genomic_to_gene(start_genomic, gene_id)
        end_gene = self.coordinate_genomic_to_gene(end_genomic - 1, gene_id)
        if transcript_model.transcript.strand == -1:
            start_gene, end_gene = end_gene, start_gene
        end_gene += 1

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
        end = self.coordinate_gene_to_genomic(feature.location.end - 1, gene_id)
        strand = gene_model.location.strand
        if strand == -1:
            start, end = end, start
        end += 1
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=strand)
        new_feature = feature.__class__(
            chrom=gene_id, location=location, attributes={}
        )
        if hasattr(feature, 'source'):
            new_feature.source = feature.source
        return new_feature

    def feature_coordinate_genomic_to_gene(self, feature:SeqFeature,
            gene_id:str) -> SeqFeature:
        """ For a given feature, converts the coordinate from genomic to gene.
        """
        gene_model = self.genes[gene_id]
        start = self.coordinate_genomic_to_gene(feature.location.start, gene_id)
        end = self.coordinate_genomic_to_gene(feature.location.end - 1, gene_id)
        strand = gene_model.location.strand
        if strand == -1:
            start, end = end, start
        end += 1
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=0)

        new_feature =  feature.__class__(
            chrom=gene_id, location=location, attributes={}
        )
        if hasattr(feature, 'source'):
            new_feature.source = feature.source

        return new_feature

    def create_gene_id_version_mapper(self) -> None:
        """ Create the a dict that keys are the unversioned gene ID, and values
        are the versioned. """
        self.gene_id_version_mapper = {}
        for versioned in self.genes.keys():
            unversioned = versioned.split('.')[0]
            if unversioned in self.gene_id_version_mapper:
                if '_PAR_Y' in versioned:
                    continue
                if '_PAR_Y' not in self.gene_id_version_mapper[unversioned]:
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
            try:
                return self.genes[gene_id]
            except KeyError as error:
                raise err.GeneNotFoundError(gene_id) from error

        if self.gene_id_version_mapper is None:
            self.create_gene_id_version_mapper()

        try:
            versioned_gene_id = self.gene_id_version_mapper[gene_id]
        except KeyError as error:
            raise err.GeneNotFoundError(gene_id) from error

        return self.genes[versioned_gene_id]

    def get_all_exons_of_gene(self, gene_id) -> List[GTFSeqFeature]:
        """ Get all exons of a gene """
        exons = set()
        gene_model = self.genes[gene_id]
        for tx_id in gene_model.transcripts:
            tx_model = self.transcripts[tx_id]
            exons.update(tx_model.exon)
        exons = list(exons)
        exons.sort()
        gene_model.exons = exons

    def find_exon_index(self, transcript_id:str, feature:GTFSeqFeature,
            coordinate:str='gene') -> int:
        """ Find the exon index of the gene.

        Args:
            gene_id (str): The gene ID.
            feature (GTFSeqFeature): The exon to look up.
            coordinate (str): The coordinate type of the given feature. Can be
                one of 'genomic' or 'gene'. Defaults to gene.
        """
        tx_model = self.transcripts[transcript_id]
        gene_id = tx_model.transcript.gene_id
        exons = tx_model.exon

        strand = tx_model.transcript.strand

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
        raise err.ExonNotFoundError(gene_id, feature)

    def find_intron_index(self, transcript_id:str, feature:GTFSeqFeature,
            coordinate:str="gene", intron_start_range:Tuple[int,int]=(0,0),
            intron_end_range:Tuple[int,int]=(0,0)) -> int:
        """ Find the intron index of the gene.

        Args:
            gene_id (str): The gene ID.
            feature (GTFSeqFeature): The intron to look up.
            coordinate (str): The coordinate type of the given feature. Can be
                one of 'genomic' or 'gene'. Defaults to gene.
        """
        tx_model = self.transcripts[transcript_id]
        exons = tx_model.exon
        strand = tx_model.transcript.strand
        gene_id = tx_model.transcript.gene_id

        intron_start_range = FeatureLocation(
            start=intron_start_range[0],
            end=intron_start_range[1] + 1
        )
        intron_end_range = FeatureLocation(
            start=intron_end_range[0],
            end=intron_end_range[1] + 1
        )

        if coordinate == 'gene':
            feature = self.feature_coordinate_gene_to_genomic(feature, gene_id)

        elif coordinate != 'genomic':
            raise ValueError("Don't know how to handle the coordinate.")

        if strand == 1:
            it = enumerate(exons)
            while True:
                i, exon = next(it, (None, None))
                if not exon:
                    break
                start_offset = feature.location.start - exon.location.end
                if start_offset in intron_start_range:
                    i, exon = next(it, (None, None))
                    if not exon:
                        break
                    end_offset = feature.location.end - exon.location.start
                    if end_offset in intron_end_range:
                        return i - 1
                    if exon.location.start >= feature.location.end:
                        return i - 1
                    break
                if exon.location.start > feature.location.end:
                    break
        elif strand == -1:
            it = reversed(list(enumerate(exons)))
            while True:
                i, exon = next(it, (None, None))
                if not exon:
                    break
                start_offset = - (feature.location.end - exon.location.start)
                if start_offset in intron_start_range:
                    i, exon = next(it, (None, None))
                    if not exon:
                        break
                    end_offset = - (feature.location.start - exon.location.end)
                    if end_offset in intron_end_range:
                        return i + 1
                    if exon.location.end <= feature.location.start:
                        return i + 1
                    break
                if exon.location.end < feature.location.start:
                    break
        raise err.IntronNotFoundError(gene_id, feature)

    def get_transcripts_with_exonic_position(self, gene_id:str, pos:int
            ) -> List[TranscriptAnnotationModel]:
        """ get all transcripts of a gene that the given genomic position is
        exonic """
        transcripts:List[TranscriptAnnotationModel] = []
        for tx_id in self.genes[gene_id].transcripts:
            tx_model = self.transcripts[tx_id]
            if tx_model.is_exonic(pos):
                transcripts.append(tx_model)
        return transcripts

    def get_transcripts_with_position(self, gene_id:str, pos:str
            ) -> List[TranscriptAnnotationModel]:
        """ Get all transcripts of a gene that contains the genomic posision in
        exon or intron """
        transcripts:List[TranscriptAnnotationModel] = []
        for tx_id in self.genes[gene_id].transcripts:
            tx_model = self.transcripts[tx_id]
            if not tx_model.exon:
                continue
            start = tx_model.transcript.location.start
            end = tx_model.transcript.location.end
            if start <= pos < end:
                transcripts.append(tx_model)
        return transcripts
