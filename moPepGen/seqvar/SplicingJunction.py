""" Splicing junction site """
from __future__ import annotations
from typing import TYPE_CHECKING, List
from moPepGen.SeqFeature import FeatureLocation
from .VariantRecord import VariantRecord


if TYPE_CHECKING:
    from moPepGen.gtf import TranscriptAnnotationModel, GenomicAnnotation
    from moPepGen.dna import DNASeqRecord

class SpliceJunction():
    """ Represents a splice junction between two genomic positions. """
    def __init__(self, upstream_start:int, upstream_end:int, downstream_start:int,
            downstream_end:int, gene_id:str, chrom:str):
        """ Constructor """
        self.upstream_start = upstream_start
        self.upstream_end = upstream_end
        self.downstream_start = downstream_start
        self.downstream_end = downstream_end
        self.gene_id = gene_id
        self.chrom = chrom

    def __str__(self) -> str:
        """ str """
        return f"<{self.__class__}: {self.upstream_start}-{self.upstream_end}:" +\
            f"{self.downstream_start}:{self.downstream_end}, "

    def align_to_transcript(self, tx_model:TranscriptAnnotationModel,
            upstream_novel:bool, downstream_novel:bool
            ) -> SpliceJunctionTranscriptAlignment:
        """ Align the junction site to a transcript.

        Args:
            - `tx_model` (TranscriptAnnotationModel): The transcript to be
              aligned to.
            - `upstream_novel` (bool): If true, the `downstream_start` must be
              matched to an exon.
            - `downstream_novel` (bool): If true, the `upstream_start` must be
              matched to an exon.
        """
        upstream_start_index = tx_model.get_exon_with_start(self.upstream_start)
        upstream_end_index = tx_model.get_exon_with_end(self.upstream_end)
        downstream_start_index = tx_model.get_exon_with_start(self.downstream_start)
        downstream_end_index = tx_model.get_exon_with_end(self.downstream_end)

        if upstream_novel and downstream_start_index == -1:
            return None

        if downstream_novel and upstream_end_index == -1:
            return None

        return SpliceJunctionTranscriptAlignment(
            junction=self,
            tx_model=tx_model,
            upstream_start_index=upstream_start_index,
            upstream_end_index=upstream_end_index,
            downstream_start_index=downstream_start_index,
            downstream_end_index=downstream_end_index,
            upstream_novel=upstream_novel,
            downstream_novel=downstream_novel
        )

    def is_novel(self, anno:GenomicAnnotation):
        """ Whether the junction is novel, i.e. not present in at least one
        transcript isoform of the gene. """
        gene_model = anno.genes[self.gene_id]
        junction = FeatureLocation(
            start=self.upstream_end, end=self.downstream_start
        )
        for tx_id in gene_model.transcripts:
            tx_model = anno.transcripts[tx_id]
            if tx_model.has_junction(junction):
                return False
        return True

class SpliceJunctionTranscriptAlignment():
    """ Alignment between a splice junction to a transcript. """
    def __init__(self, junction:SpliceJunction, tx_model:TranscriptAnnotationModel,
            upstream_start_index:int, upstream_end_index:int,
            downstream_start_index:int, downstream_end_index:int,
            upstream_novel:bool, downstream_novel:bool):
        """ Constructor """
        self.junction = junction
        self.tx_model = tx_model
        self.upstream_start_index = upstream_start_index
        self.upstream_end_index = upstream_end_index
        self.downstream_start_index = downstream_start_index
        self.downstream_end_index = downstream_end_index
        self.upstream_novel = upstream_novel
        self.downstream_novel = downstream_novel

    def get_interjacent_exons(self) -> List[int]:
        """ Get interjacent exons between upstream end and downstream start. """
        if not self.upstream_end_index and not self.downstream_start_index:
            raise ValueError(
                "Neither upstream end nor downstream start has a matched exon."
            )

        interjacent = []
        if self.upstream_end_index > -1:
            if self.upstream_end_index + 1 == len(self.tx_model.exon):
                return interjacent
            iter_exon = range(self.upstream_end_index + 1, len(self.tx_model.exon))
            is_reversed = False
        else:
            if self.downstream_start_index == 0:
                return interjacent
            iter_exon = range(self.downstream_start_index - 1, -1, -1)
            is_reversed = True

        for i in iter_exon:
            exon = self.tx_model.exon[i]
            if self.junction.upstream_end <= exon.location.start \
                    < exon.location.end <= self.junction.downstream_start:
                interjacent.append(i)
            if exon.location.end <= self.junction.upstream_end \
                    or exon.location.start >= self.junction.downstream_start:
                break

        if is_reversed:
            interjacent = list(reversed(interjacent))
        return interjacent

    def get_upstream_end_spanning(self) -> int:
        """ Get exons that are spanning over the upstream end position. """
        if self.downstream_start_index == -1:
            return self.tx_model.get_exon_containing(self.junction.upstream_end - 1)

        i = self.downstream_start_index - 1
        while i >= 0:
            exon = self.tx_model.exon[i]
            if self.junction.upstream_end - 1 in exon.location:
                return i
            i -= 1
        return -1

    def get_downstream_start_spanning(self) -> int:
        """ Get exons that are spanning over the downstream start position. """
        if self.upstream_end_index == -1:
            return self.tx_model.get_exon_containing(self.junction.downstream_start)

        i = self.upstream_end_index + 1
        while i < len(self.tx_model.exon):
            exon = self.tx_model.exon[i]
            if self.junction.downstream_start in exon.location:
                return i
            i += 1
        return -1

    def create_upstream_deletion(self, spanning:int, interjacent:List[int],
            anno:GenomicAnnotation, gene_seq:DNASeqRecord, var_id:str
            ) -> VariantRecord:
        """ Create an upstream deletion variant. """
        gene_id = self.junction.gene_id
        gene_model = anno.genes[gene_id]
        tx_id = self.tx_model.transcript.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = self.tx_model

        if tx_model.exon[spanning].location.end == self.junction.upstream_end:
            genomic_start = tx_model.exon[interjacent[0]].location.start
        else:
            genomic_start = self.junction.upstream_end
        start = anno.coordinate_genomic_to_gene(genomic_start, gene_id)

        if interjacent:
            genomic_end = tx_model.exon[interjacent[-1]].location.end
        else:
            genomic_end = tx_model.exon[spanning].location.end
        end = anno.coordinate_genomic_to_gene(genomic_end - 1, gene_id)

        if strand == -1:
            start, end = end, start
        end += 1

        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

        location = FeatureLocation(seqname=gene_id, start=start, end=end)
        alt = '<DEL>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Deletion'
        return VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_downstream_deletion(self, spanning:int, interjacent:List[int],
            anno:GenomicAnnotation, gene_seq:DNASeqRecord, var_id:str
            ) -> VariantRecord:
        """ Create a downstream deletion variant. """
        gene_id = self.junction.gene_id
        gene_model = anno.genes[gene_id]
        tx_id = self.tx_model.transcript.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = self.tx_model

        if interjacent:
            genomic_start = tx_model.exon[interjacent[0]].location.start
        else:
            genomic_start = tx_model.exon[spanning].location.start
        start = anno.coordinate_genomic_to_gene(genomic_start, gene_id)

        if self.junction.downstream_start == tx_model.exon[spanning].location.start:
            genomic_end = tx_model.exon[interjacent[-1]].location.end
        else:
            genomic_end = self.junction.downstream_start
        end = anno.coordinate_genomic_to_gene(genomic_end - 1, gene_id)

        if strand == -1:
            start, end = end, start
        end += 1

        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

        location = FeatureLocation(seqname=gene_id, start=start, end=end)
        alt = '<DEL>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Deletion'
        return VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_upstream_substitution(self, interjacent:List[int],
            anno:GenomicAnnotation, gene_seq:DNASeqRecord, var_id:str
            ) -> VariantRecord:
        """ Create an upstream substitution variant. """
        gene_id = self.junction.gene_id
        gene_model = anno.genes[gene_id]
        tx_id = self.tx_model.transcript.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = self.tx_model

        genomic_start = tx_model.exon[interjacent[0]].location.start
        start = anno.coordinate_genomic_to_gene(genomic_start, gene_id)
        genomic_end = tx_model.exon[interjacent[-1]].location.end
        end = anno.coordinate_genomic_to_gene(genomic_end - 1, gene_id)

        if interjacent[0] > 0:
            genomic_donor_start = tx_model.exon[interjacent[0] - 1].location.end
            if self.junction.upstream_start is not None:
                genomic_donor_start = max(genomic_donor_start, self.junction.upstream_start)
        else:
            if self.junction.upstream_start is None:
                raise ValueError('`upstream_start` is missing.')
            genomic_donor_start = self.junction.upstream_start
        donor_start = anno.coordinate_genomic_to_gene(genomic_donor_start, gene_id)
        genomic_donor_end = self.junction.upstream_end
        donor_end = anno.coordinate_genomic_to_gene(genomic_donor_end - 1, gene_id)

        if strand == -1:
            start, end = end, start
            donor_start, donor_end = donor_end, donor_start

        end += 1
        donor_end += 1

        ref = str(gene_seq.seq[start])
        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

        location = FeatureLocation(seqname=gene_id, start=start, end=end)
        alt = '<SUB>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'DONOR_GENE_ID': gene_id,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Substitution'
        return VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_downstream_substitution(self, interjacent:List[int],
            anno:GenomicAnnotation, gene_seq:DNASeqRecord, var_id:str
            ) -> VariantRecord:
        """ Create a downstream substitution variant. """
        gene_id = self.junction.gene_id
        gene_model = anno.genes[gene_id]
        tx_id = self.tx_model.transcript.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = self.tx_model

        genomic_start = tx_model.exon[interjacent[0]].location.start
        start = anno.coordinate_genomic_to_gene(genomic_start, gene_id)
        genomic_end = tx_model.exon[interjacent[-1]].location.end
        end = anno.coordinate_genomic_to_gene(genomic_end - 1, gene_id)

        genomic_donor_start = self.junction.downstream_start
        donor_start = anno.coordinate_genomic_to_gene(genomic_donor_start, gene_id)

        if interjacent[-1] < len(tx_model.exon) - 1:
            genomic_donor_end = tx_model.exon[interjacent[-1] + 1].location.start
            if self.junction.downstream_end is not None:
                genomic_donor_end = min(genomic_donor_end, self.junction.downstream_end)
        else:
            if self.junction.downstream_end is None:
                raise ValueError('`downstream_start` is missing.')
            genomic_donor_end = self.junction.downstream_end
        donor_end = anno.coordinate_genomic_to_gene(genomic_donor_end - 1, gene_id)

        if strand == -1:
            start, end = end, start
            donor_start, donor_end = donor_end, donor_start

        end += 1
        donor_end += 1

        ref = str(gene_seq.seq[start])
        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

        location = FeatureLocation(seqname=gene_id, start=start, end=end)
        alt = '<SUB>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'DONOR_GENE_ID': gene_id,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Substitution'
        return VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_upstream_insertion(self, anno:GenomicAnnotation,
            gene_seq:DNASeqRecord, var_id:str) -> VariantRecord:
        """ Create an upstream insertion variant. """
        gene_id = self.junction.gene_id
        gene_model = anno.genes[gene_id]
        tx_id = self.tx_model.transcript.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = self.tx_model

        if self.downstream_start_index <= 0:
            raise ValueError("`downstream_start_index` should not be > 0.")

        if strand == 1:
            insert_pos_genomic = tx_model.exon[self.downstream_start_index - 1].location.end - 1
            donor_start = tx_model.exon[self.downstream_start_index - 1].location.end
            if self.junction.downstream_end is not None:
                donor_start = max(donor_start, self.junction.upstream_start)
            donor_end = self.junction.upstream_end - 1
        else:
            insert_pos_genomic = self.junction.downstream_start
            donor_start = self.junction.upstream_end - 1
            donor_end = tx_model.exon[self.downstream_start_index - 1].location.end
            if self.junction.upstream_start is not None:
                donor_end = max(donor_end, self.junction.upstream_start)

        insert_pos = anno.coordinate_genomic_to_gene(insert_pos_genomic, gene_id)

        donor_start = anno.coordinate_genomic_to_gene(donor_start, gene_id)
        donor_end = anno.coordinate_genomic_to_gene(donor_end, gene_id) + 1

        ref = str(gene_seq.seq[insert_pos])
        genomic_position = f'{chrom}:{insert_pos_genomic + 1}:{insert_pos_genomic + 2}'

        location = FeatureLocation(
            seqname=gene_id,
            start=insert_pos,
            end=insert_pos + 1
        )
        alt = '<Ins>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'DONOR_GENE_ID': gene_id,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Insertion'
        return VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_downstream_insertion(self, anno:GenomicAnnotation,
            gene_seq:DNASeqRecord, var_id:str) -> VariantRecord:
        """ Create a downstream insertion variant. """
        gene_id = self.junction.gene_id
        gene_model = anno.genes[gene_id]
        tx_id = self.tx_model.transcript.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = self.tx_model

        if self.upstream_end_index == -1:
            raise ValueError("`upstream_start_index` should not be -1.")

        if self.upstream_start_index == len(tx_model.exon) - 1:
            raise ValueError("`downstream_start_index` should not be the last one.")

        if strand == 1:
            insert_pos_genomic = self.junction.upstream_end - 1
            donor_start = self.junction.downstream_start
            donor_end = tx_model.exon[self.upstream_end_index + 1].location.start - 1
            if self.junction.downstream_end is not None:
                donor_end = min(donor_end, self.junction.downstream_end - 1)
        else:
            insert_pos_genomic = tx_model.exon[self.upstream_end_index + 1].location.start
            donor_start = tx_model.exon[self.upstream_end_index + 1].location.start - 1
            if self.junction.downstream_end is not None:
                donor_start = min(donor_start, self.junction.downstream_end - 1)
            donor_end = self.junction.downstream_start
        insert_pos = anno.coordinate_genomic_to_gene(insert_pos_genomic, gene_id)

        donor_start = anno.coordinate_genomic_to_gene(donor_start, gene_id)
        donor_end = anno.coordinate_genomic_to_gene(donor_end, gene_id) + 1

        ref = str(gene_seq.seq[insert_pos])
        genomic_position = f'{chrom}:{insert_pos_genomic + 1}:{insert_pos_genomic + 2}'

        location = FeatureLocation(
            seqname=gene_id,
            start=insert_pos,
            end=insert_pos + 1
        )
        alt = '<Ins>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'DONOR_GENE_ID': gene_id,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Insertion'
        return VariantRecord(location, ref, alt, _type, var_id, attrs)

    def convert_to_variant_records(self, anno:GenomicAnnotation,
            gene_seq:DNASeqRecord, var_id:str):
        """ Convert splicing to variant records. """
        variants = []
        interjacent = self.get_interjacent_exons()

        if self.upstream_novel:
            if self.upstream_end_index == -1 or len(interjacent) > 0:
                spanning = self.get_upstream_end_spanning()
                if spanning > -1:
                    v = self.create_upstream_deletion(
                        spanning, interjacent, anno, gene_seq, var_id
                    )
                    variants.append(v)
                elif len(interjacent) > 0:
                    v = self.create_upstream_substitution(
                        interjacent, anno, gene_seq, var_id
                    )
                    variants.append(v)
                elif self.downstream_start_index > 0:
                    v = self.create_upstream_insertion(
                        anno, gene_seq, var_id
                    )
                    variants.append(v)

        if self.downstream_novel:
            if self.downstream_start_index == -1 or len(interjacent) > 0:
                spanning = self.get_downstream_start_spanning()
                if spanning > -1:
                    v = self.create_downstream_deletion(
                        spanning, interjacent, anno, gene_seq, var_id
                    )
                    variants.append(v)
                elif len(interjacent) > 0:
                    v = self.create_downstream_substitution(
                        interjacent, anno, gene_seq, var_id
                    )
                    variants.append(v)
                elif -1 < self.downstream_end_index < len(self.tx_model.exon) - 1:
                    v = self.create_downstream_insertion(
                        anno, gene_seq, var_id
                    )
                    variants.append(v)

        if not self.downstream_novel and not self.upstream_novel:
            strand = self.tx_model.transcript.strand
            is_upstream_aligned = (strand == 1 and self.upstream_end_index != -1) \
                    or (strand == -1 and self.downstream_start_index != -1)
            tx_start = self.tx_model.transcript.location.start
            tx_end = self.tx_model.transcript.location.end
            is_after_tx_start = \
                (strand == 1 and tx_start < self.junction.upstream_end) \
                or (strand == -1 and tx_end > self.junction.downstream_start + 1)
            if is_upstream_aligned and is_after_tx_start:
                if strand == 1:
                    if self.downstream_start_index == -1 or len(interjacent) > 0:
                        spanning = self.get_downstream_start_spanning()
                        if spanning > -1:
                            v = self.create_downstream_deletion(
                                spanning, interjacent, anno, gene_seq, var_id
                            )
                            variants.append(v)
                        elif len(interjacent) > 0:
                            v = self.create_downstream_substitution(
                                interjacent, anno, gene_seq, var_id
                            )
                            variants.append(v)
                else:
                    spanning = self.get_upstream_end_spanning()
                    if self.upstream_end_index == -1 or len(interjacent) > 0:
                        if spanning > -1:
                            v = self.create_upstream_deletion(
                                spanning, interjacent, anno, gene_seq, var_id
                            )
                            variants.append(v)
                        elif len(interjacent) > 0:
                            v = self.create_upstream_substitution(
                                interjacent, anno, gene_seq, var_id
                            )
                            variants.append(v)

        return variants
