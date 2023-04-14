""" Module for rMATS parser """
from __future__ import annotations
from typing import TYPE_CHECKING, List, Dict
from moPepGen import gtf, dna, seqvar
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord


if TYPE_CHECKING:
    from moPepGen.gtf import TranscriptAnnotationModel
    from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature

class SETranscriptAlignment():
    """ Alignment of an SE coordinates to an transcript exons. """
    def __init__(self, se_record:SERecord, tx_model:TranscriptAnnotationModel,
            upstream_index:int, downstream_index:int, target_indices:List[int]):
        """ Constructor """
        self.se_record = se_record
        self.tx_model = tx_model
        self.upstream_index = upstream_index
        self.downstream_index = downstream_index
        self.target_indices = target_indices

    def get_upstream_interjacent_exons(self) -> List[GTFSeqFeature]:
        """ Get upstream interjacent exons between the annotated upstream exon
        and the exon being skipped. Returned values are exon indices. """
        interjacent_exons = []
        for exon in self.tx_model.exon[self.upstream_index + 1:]:
            if exon.location.start < self.se_record.exon_start:
                interjacent_exons.append(exon)
            if exon.location.start >= self.se_record.exon_start:
                break
        return interjacent_exons

    def get_downstream_interjacent_exons(self) -> List[GTFSeqFeature]:
        """ Get downstream interjacent exons between the annotated upstream exon
        and the exon being skipped. Returned values are exon indices. """
        interjacent_exons = []
        for exon in reversed(self.tx_model.exon[:self.downstream_index]):
            if exon.location.end > self.se_record.exon_end:
                interjacent_exons.append(exon)
            if exon.location.end <= self.se_record.exon_end:
                break
        return list(reversed(interjacent_exons))

class SERecord(RMATSRecord):
    """ Skipped Exon """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str, exon_start:int,
            exon_end:int, upstream_exon_start:int, upstream_exon_end:int,
            downstream_exon_start:int, downstream_exon_end:int,
            ijc_sample_1:int, sjc_sample_1:int, ijc_sample_2:int,
            sjc_sample_2:int, inc_form_len:int, skip_form_len:int,
            pvalue:float, fdr:float):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.upstream_exon_start = upstream_exon_start
        self.upstream_exon_end = upstream_exon_end
        self.downstream_exon_start = downstream_exon_start
        self.downstream_exon_end = downstream_exon_end
        self.ijc_sample_1 = ijc_sample_1
        self.sjc_sample_1 = sjc_sample_1
        self.ijc_sample_2 = ijc_sample_2
        self.sjc_sample_2 = sjc_sample_2
        self.inc_form_len = inc_form_len
        self.skip_form_len = skip_form_len
        self.pvalue = pvalue
        self.fdr = fdr

    @classmethod
    def readline(cls, line:str) -> SERecord:
        """ Create a SERecord from a line of the rMATS output """
        fields = line.rstrip().split('\t')
        return cls(
            gene_id=fields[1].strip('"'),
            gene_symbol=fields[2].strip('"'),
            chrom=fields[3],
            exon_start=int(fields[5]),
            exon_end=int(fields[6]),
            upstream_exon_start=int(fields[7]),
            upstream_exon_end=int(fields[8]),
            downstream_exon_start=int(fields[9]),
            downstream_exon_end=int(fields[10]),
            ijc_sample_1=int(fields[12]),
            sjc_sample_1=int(fields[13]),
            ijc_sample_2=None if fields[14] == '' else int(fields[14]),
            sjc_sample_2=None if fields[15] == '' else int(fields[15]),
            inc_form_len=int(fields[16]),
            skip_form_len=int(fields[17]),
            pvalue=None if fields[18] == 'NA' else float(fields[18]),
            fdr=None if fields[19] == 'NA' else float(fields[19])
        )

    def align_to_transcript(self, tx_model:gtf.TranscriptAnnotationModel,
            se_loc:FeatureLocation) -> SETranscriptAlignment:
        """ Align the SE record to a transcript """
        upstream_index = tx_model.get_exon_with_end(self.upstream_exon_end)
        if upstream_index == -1:
            return None
        downstream_index = tx_model.get_exon_with_start(
            self.downstream_exon_start, upstream_index + 1
        )
        if downstream_index == -1:
            return None

        target_exons = tx_model.get_exon_inner(se_loc, upstream_index)

        return SETranscriptAlignment(
            se_record=self,
            tx_model=tx_model,
            upstream_index=upstream_index,
            downstream_index=downstream_index,
            target_indices=target_exons
        )

    def create_variant_id(self, anno:gtf.GenomicAnnotation) -> str:
        """ Create variant ID """
        uee = anno.coordinate_genomic_to_gene(self.upstream_exon_end, self.gene_id)
        es = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        ee = anno.coordinate_genomic_to_gene(self.exon_end - 1, self.gene_id)
        des = anno.coordinate_genomic_to_gene(self.downstream_exon_start - 1, self.gene_id)

        gene_model = anno.genes[self.gene_id]
        if gene_model.strand == -1:
            es, ee = ee, es
            uee, des = des, uee
        ee += 1
        des += 1

        return f"SE_{uee}-{es}-{ee}-{des}"

    def has_novel_splicing(self, anno:gtf.GenomicAnnotation) -> bool:
        """ Checks if the AS event is novel, i.e. whether both retained and
        skipped version already exist in reference. """
        gene_model = anno.genes[self.gene_id]

        upstream_junction = FeatureLocation(
            start=self.upstream_exon_end, end=self.exon_start
        )
        downstream_junction = FeatureLocation(
            start=self.exon_end, end=self.downstream_exon_start
        )
        skip_junction = FeatureLocation(
            start=self.upstream_exon_end, end=self.downstream_exon_start
        )

        has_retained = False
        has_skipped = False

        for tx_id in gene_model.transcripts:
            tx_model = anno.transcripts[tx_id]
            if tx_model.has_junction(upstream_junction) \
                    and tx_model.has_junction(downstream_junction):
                has_retained = True
            if tx_model.has_junction(skip_junction):
                has_skipped = True

        return not has_retained or not has_skipped

    def create_target_deletion(self, alignment:SETranscriptAlignment,
            anno:gtf.GenomicAnnotation, gene_seq:dna.DNASeqRecord, var_id:str
            ) -> seqvar.VariantRecord:
        """ Create a deletion of the target exons """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        start = alignment.tx_model.exon[alignment.upstream_index + 1].location.start
        start = anno.coordinate_genomic_to_gene(start, self.gene_id)
        end = alignment.tx_model.exon[alignment.downstream_index - 1].location.end
        end = anno.coordinate_genomic_to_gene(end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        alt = '<DEL>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Deletion'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_upstream_interjacent_deletion(self, alignment:SETranscriptAlignment,
            upstream_interjacent:List[GTFSeqFeature], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ Create a deletion of the upstream interjacent exons. """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        start = upstream_interjacent[0].location.start
        start = anno.coordinate_genomic_to_gene(start, self.gene_id)
        end = self.exon_start
        end = anno.coordinate_genomic_to_gene(end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        alt = '<DEL>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Deletion'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_upstream_interjacent_substitution(self, alignment:SETranscriptAlignment,
            upstream_interjacent:List[GTFSeqFeature], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ Create a substitution of the upstream interjacent exons with the
        skipped exon when it is included. """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        start = upstream_interjacent[0].location.start
        start = anno.coordinate_genomic_to_gene(start, self.gene_id)
        end = self.exon_start
        end = anno.coordinate_genomic_to_gene(end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        donor_start = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        donor_end = anno.coordinate_genomic_to_gene(self.exon_end - 1, self.gene_id)
        if gene_model.strand == -1:
            donor_start, donor_end = donor_end, donor_start
        donor_end += 1

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        alt = '<SUB>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'DONOR_GENE_ID': self.gene_id,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Substitution'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_downstream_interjacent_deletion(self, alignment:SETranscriptAlignment,
            downstream_interjacent:List[GTFSeqFeature], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ Create a deletion of the downstream interjacent exons. """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        start = self.exon_end
        start = anno.coordinate_genomic_to_gene(start, self.gene_id)
        end = downstream_interjacent[-1].location.end
        end = anno.coordinate_genomic_to_gene(end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        alt = '<DEL>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Deletion'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_downstream_interjacent_substitution(self, alignment:SETranscriptAlignment,
            downstream_interjacent:List[GTFSeqFeature], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ Create a substitution of the downstream interjacent exons with the
        target skipped exon when it is included. """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        start = self.exon_end
        start = anno.coordinate_genomic_to_gene(start, self.gene_id)
        end = downstream_interjacent[-1].location.end
        end = anno.coordinate_genomic_to_gene(end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        donor_start = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        donor_end = anno.coordinate_genomic_to_gene(self.exon_end - 1, self.gene_id)
        if gene_model.strand == -1:
            donor_start, donor_end = donor_end, donor_start
        donor_end += 1

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        alt = '<SUB>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'DONOR_GENE_ID': self.gene_id,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Substitution'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_target_insertion(self, alignment:SETranscriptAlignment,
            anno:gtf.GenomicAnnotation, gene_seq:dna.DNASeqRecord, var_id:str
            ) -> seqvar.VariantRecord:
        """ Create a insertion of the target skipped exon when it is included. """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom

        if gene_model.strand == 1:
            insert_position = anno.coordinate_genomic_to_gene(
                index=self.upstream_exon_end - 1,
                gene=self.gene_id
            )
        else:
            insert_position = anno.coordinate_genomic_to_gene(
                index=self.downstream_exon_start,
                gene=self.gene_id
            )

        donor_start = self.exon_start
        donor_start = anno.coordinate_genomic_to_gene(donor_start, self.gene_id)
        donor_end = self.exon_end
        donor_end = anno.coordinate_genomic_to_gene(donor_end - 1, self.gene_id)
        if gene_model.strand == -1:
            donor_start, donor_end = donor_end, donor_start
        donor_end += 1

        ref = str(gene_seq.seq[insert_position])

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(
            seqname=self.gene_id,
            start=insert_position,
            end=insert_position + 1
        )
        alt = '<Ins>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'DONOR_GENE_ID': self.gene_id,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Insertion'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_target_substitution(self, alignment:SETranscriptAlignment,
            anno:gtf.GenomicAnnotation, gene_seq:dna.DNASeqRecord, var_id:str
            ) -> seqvar.VariantRecord:
        """ Create a substitution of the upstream and/or downstream interjacent
        exons with the skipped exon when it is included. """
        interjacent_exons = alignment.tx_model.exon[
            alignment.upstream_index + 1:alignment.downstream_index]
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        start = interjacent_exons[0].location.start
        start = anno.coordinate_genomic_to_gene(start, self.gene_id)
        end = interjacent_exons[-1].location.end
        end = anno.coordinate_genomic_to_gene(end - 1, self.gene_id)
        if gene_model.strand == -1:
            start, end = end, start
        end += 1
        ref = str(gene_seq.seq[start])

        donor_start = anno.coordinate_genomic_to_gene(self.exon_start, self.gene_id)
        donor_end = anno.coordinate_genomic_to_gene(self.exon_end - 1, self.gene_id)
        if gene_model.strand == -1:
            donor_start, donor_end = donor_end, donor_start
        donor_end += 1

        genomic_position = f'{chrom}:{self.exon_start+1}:{self.exon_end}'

        location = FeatureLocation(seqname=self.gene_id, start=start, end=end)
        alt = '<SUB>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': start,
            'END': end,
            'DONOR_START': donor_start,
            'DONOR_END': donor_end,
            'DONOR_GENE_ID': self.gene_id,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Substitution'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        """ Convert to variant records """
        variants = []
        if not self.has_novel_splicing(anno):
            return variants

        gene_model = anno.genes[self.gene_id]
        transcript_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])

        se_loc = FeatureLocation(self.exon_start, self.exon_end)

        tx_alignments:Dict[str, SETranscriptAlignment] = {}

        for tx_id in transcript_ids:
            model = anno.transcripts[tx_id]

            alignment = self.align_to_transcript(model, se_loc)
            if alignment:
                tx_alignments[tx_id] = alignment

        var_id = self.create_variant_id(anno)

        # For SE, the skipped is 'skipped', and inclusion is 'inclusion'
        # what a nonsense comment
        for tx_id, alignment in tx_alignments.items():
            if alignment.target_indices and self.sjc_sample_1 >= min_sjc:
                record = self.create_target_deletion(alignment, anno, gene_seq, var_id)
                variants.append(record)

            upstream_interjacent = alignment.get_upstream_interjacent_exons()
            if upstream_interjacent:
                if alignment.target_indices:
                    # This is a Deletion
                    if self.ijc_sample_1 >= min_ijc:
                        record = self.create_upstream_interjacent_deletion(
                            alignment, upstream_interjacent, anno, gene_seq,
                            var_id
                        )
                        variants.append(record)
                else:
                    # This is a Substitution
                    if self.sjc_sample_1 >= min_sjc:
                        record = self.create_upstream_interjacent_substitution(
                            alignment, upstream_interjacent, anno, gene_seq,
                            var_id
                        )
                        variants.append(record)

            downstream_interjacent = alignment.get_downstream_interjacent_exons()
            if downstream_interjacent:
                if alignment.target_indices:
                    if self.ijc_sample_1 >= min_ijc:
                        record = self.create_downstream_interjacent_deletion(
                            alignment, downstream_interjacent, anno, gene_seq,
                            var_id
                        )
                        variants.append(record)
                else:
                    if self.sjc_sample_1 >= min_sjc:
                        record = self.create_downstream_interjacent_substitution(
                            alignment, downstream_interjacent, anno, gene_seq,
                            var_id
                        )
                        variants.append(record)

            if not alignment.target_indices and self.ijc_sample_1 >= min_ijc:
                if alignment.upstream_index + 1 == alignment.downstream_index:
                    # No interjacent exons
                    record = self.create_target_insertion(alignment, anno, gene_seq, var_id)
                else:
                    record = self.create_target_substitution(alignment, anno, gene_seq, var_id)
                variants.append(record)
        return variants
