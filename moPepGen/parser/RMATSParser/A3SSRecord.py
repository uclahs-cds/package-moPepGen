""" Module for rMATS parser """
from __future__ import annotations
from typing import TYPE_CHECKING, List, Tuple, Dict
from moPepGen import gtf, dna, seqvar
from moPepGen.SeqFeature import FeatureLocation
from .RMATSRecord import RMATSRecord


if TYPE_CHECKING:
    from moPepGen.gtf import TranscriptAnnotationModel
    from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature

class A3SSTranscriptAlignment():
    """ Alignment of an SE coordinates to an transcript exons. """
    def __init__(self, as_record:A3SSRecord, tx_model:TranscriptAnnotationModel,
            upstream_index:int, long_index:int, short_index:int):
        """ Constructor """
        self.as_record = as_record
        self.tx_model = tx_model
        self.upstream_index = upstream_index
        self.long_index = long_index
        self.short_index = short_index

    def get_long_interjacent_exons(self) -> List[int]:
        """ Get the index of interjacent exons between upstream exon and the
        long exon. """
        interjacent = []
        strand = self.tx_model.transcript.strand
        if strand == 1:
            for i in range(self.upstream_index + 1, len(self.tx_model.exon)):
                exon = self.tx_model.exon[i]
                if exon.location.start < self.as_record.long_exon_start:
                    interjacent.append(i)
                if exon.location.start >= self.as_record.long_exon_start:
                    break
        else:
            for i in reversed(range(0, self.upstream_index)):
                exon = self.tx_model.exon[i]
                if exon.location.end > self.as_record.long_exon_end:
                    interjacent.append(i)
                if exon.location.end <= self.as_record.long_exon_end:
                    break
            interjacent = list(reversed(interjacent))
        return interjacent

    def get_short_interjacent_exons(self) -> List[GTFSeqFeature]:
        """ Get the index of interjacent exons between upstream exon and the
        short exon. """
        interjacent = []
        strand = self.tx_model.transcript.strand
        if strand == 1:
            for i in range(self.upstream_index + 1, len(self.tx_model.exon)):
                exon = self.tx_model.exon[i]
                if exon.location.start < self.as_record.short_exon_start:
                    interjacent.append(i)
                if exon.location.start >= self.as_record.short_exon_start:
                    break
        else:
            for i in reversed(range(0, self.upstream_index)):
                exon = self.tx_model.exon[i]
                if exon.location.end > self.as_record.short_exon_end:
                    interjacent.append(i)
                if exon.location.end <= self.as_record.short_exon_end:
                    break
            interjacent = list(reversed(interjacent))
        return interjacent

    def get_long_spanning_exon(self) -> int:
        """ """
        strand = self.tx_model.transcript.strand
        if strand == 1:
            for i in range(self.upstream_index + 1, len(self.tx_model.exon)):
                exon = self.tx_model.exon[i]
                if self.as_record.long_exon_start in exon.location:
                    return i
        else:
            for i in reversed(range(0, self.upstream_index)):
                exon = self.tx_model.exon[i]
                if self.as_record.long_exon_start in exon.location:
                    return i
        return -1

    def get_short_spanning_exon(self) -> int:
        """ """
        strand = self.tx_model.transcript.strand
        if strand == 1:
            for i in range(self.upstream_index + 1, len(self.tx_model.exon)):
                exon = self.tx_model.exon[i]
                if self.as_record.short_exon_start in exon.location:
                    return i
        else:
            for i in reversed(range(0, self.upstream_index)):
                exon = self.tx_model.exon[i]
                if self.as_record.short_exon_start in exon.location:
                    return i
        return -1

class A3SSRecord(RMATSRecord):
    """ Alternative 3' Splicing Site """
    def __init__(self, gene_id:str, gene_symbol:str, chrom:str,
            long_exon_start:int, long_exon_end:int, short_exon_start:int,
            short_exon_end:int, flanking_exon_start, flanking_exon_end:int,
            ijc_sample_1:int, sjc_sample_1:int, ijc_sample_2:int,
            sjc_sample_2:int, inc_form_len:int, skip_form_len:int,
            pvalue:float, fdr:float):
        """ Constructor """
        super().__init__(gene_id, gene_symbol, chrom)
        self.long_exon_start = long_exon_start
        self.long_exon_end = long_exon_end
        self.short_exon_start = short_exon_start
        self.short_exon_end = short_exon_end
        self.flanking_exon_start = flanking_exon_start
        self.flanking_exon_end = flanking_exon_end
        self.ijc_sample_1 = ijc_sample_1
        self.sjc_sample_1 = sjc_sample_1
        self.ijc_sample_2 = ijc_sample_2
        self.sjc_sample_2 = sjc_sample_2
        self.inc_form_len = inc_form_len
        self.skip_form_len = skip_form_len
        self.pvalue = pvalue
        self.fdr = fdr

    @classmethod
    def readline(cls, line:str) -> A3SSRecord:
        """ Create a A5SSRecord from a line of the rMASTS output """
        fields = line.rstrip().split('\t')
        return cls(
            gene_id=fields[1].strip('"'),
            gene_symbol=fields[2].strip('"'),
            chrom=fields[3],
            long_exon_start=int(fields[5]),
            long_exon_end=int(fields[6]),
            short_exon_start=int(fields[7]),
            short_exon_end=int(fields[8]),
            flanking_exon_start=int(fields[9]),
            flanking_exon_end=int(fields[10]),
            ijc_sample_1=int(fields[12]),
            sjc_sample_1=int(fields[13]),
            ijc_sample_2=None if fields[14] == '' else int(fields[14]),
            sjc_sample_2=None if fields[15] == '' else int(fields[15]),
            inc_form_len=int(fields[16]),
            skip_form_len=int(fields[17]),
            pvalue=None if fields[18] == 'NA' else float(fields[18]),
            fdr=None if fields[19] == 'NA' else float(fields[19])
        )

    def align_to_transcript(self, tx_model:gtf.TranscriptAnnotationModel
            ) -> A3SSTranscriptAlignment:
        """ Align the A3SS to a transcript """
        strand = tx_model.transcript.strand
        if strand == 1:
            upstream_index = tx_model.get_exon_with_end(self.flanking_exon_end)
        else:
            upstream_index = tx_model.get_exon_with_start(self.flanking_exon_start)

        if upstream_index == -1:
            return None

        if strand == 1:
            long_index = tx_model.get_exon_with_start(self.long_exon_start)
            short_index = tx_model.get_exon_with_start(self.short_exon_start)
        else:
            long_index = tx_model.get_exon_with_end(self.long_exon_end)
            short_index = tx_model.get_exon_with_end(self.short_exon_end)

        return A3SSTranscriptAlignment(
            as_record=self,
            tx_model=tx_model,
            upstream_index=upstream_index,
            long_index=long_index,
            short_index=short_index
        )

    def create_variant_id(self, anno:gtf.GenomicAnnotation) -> str:
        """ Create variant ID """
        if anno.genes[self.gene_id].strand == 1:
            uee = anno.coordinate_genomic_to_gene(self.flanking_exon_end, self.gene_id)
            les = anno.coordinate_genomic_to_gene(self.long_exon_start, self.gene_id)
            ses = anno.coordinate_genomic_to_gene(self.short_exon_start, self.gene_id)
        else:
            uee = anno.coordinate_genomic_to_gene(self.flanking_exon_start, self.gene_id)
            les = anno.coordinate_genomic_to_gene(self.long_exon_end, self.gene_id)
            ses = anno.coordinate_genomic_to_gene(self.short_exon_end, self.gene_id)

        return f"A5SS_{uee}-{les}-{ses}"

    def has_novel_splicing(self, anno:gtf.GenomicAnnotation) -> bool:
        """ Checks if the AS event is novel, i.e. whether both short and long
        version exist in reference already. """
        gene_model = anno.genes[self.gene_id]

        if gene_model.strand == 1:
            short_junction = FeatureLocation(
                start=self.flanking_exon_end, end=self.short_exon_start
            )
            long_junction = FeatureLocation(
                start=self.flanking_exon_end, end=self.long_exon_start
            )
        else:
            short_junction = FeatureLocation(
                start=self.short_exon_end, end=self.flanking_exon_start
            )
            long_junction = FeatureLocation(
                start=self.long_exon_end, end=self.flanking_exon_start
            )

        has_short = False
        has_long = False

        for tx_id in gene_model.transcripts:
            tx_model = anno.transcripts[tx_id]
            if tx_model.has_junction(short_junction):
                has_short = True
            if tx_model.has_junction(long_junction):
                has_long = True

        return not has_short or not has_long

    def create_long_deletion(self, alignment:A3SSTranscriptAlignment,
            spanning:int, interjacent:List[int], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = alignment.tx_model

        if strand == 1:
            if interjacent:
                genomic_start = tx_model.exon[interjacent[0]].location.start
            else:
                genomic_start = tx_model.exon[spanning].location.start
            start = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id)
            genomic_end = min(self.long_exon_start, tx_model.exon[spanning].location.end)
            end = anno.coordinate_gene_to_genomic(genomic_end - 1, self.gene_id) + 1
        else:
            if interjacent:
                genomic_end = tx_model.exon[interjacent[-1]].location.end
            else:
                genomic_end = tx_model.exon[spanning].location.end
            start = anno.coordinate_gene_to_genomic(genomic_end - 1, self.gene_id)
            genomic_start = max(self.long_exon.end, tx_model.exon[spanning].location.start)
            end = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id) + 1
        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

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

    def create_long_substitution(self, alignment:A3SSTranscriptAlignment,
            interjacent:List[int], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = alignment.tx_model

        if strand == 1:
            genomic_start = tx_model.exon[interjacent[0]].location.start
            start = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id)
            genomic_end = tx_model.exon[interjacent[-1]].location.end
            end = anno.coordinate_genomic_to_gene(genomic_end - 1, self.gene_id) + 1

            genomic_donor_start = self.long_exon_start
            donor_start = anno.coordinate_genomic_to_gene(genomic_donor_start, self.gene_id)
            # if the transcript does not have a downstream exon, the entire
            # long exon will be inserted.
            to_insert_full_exon = interjacent[-1] + 1 == len(tx_model.exon) \
                or tx_model.exon[interjacent[-1] + 1].location.start >= self.long_exon_end

            if to_insert_full_exon:
                genomic_donor_end = self.long_exon_end
            else:
                genomic_donor_end = tx_model.exon[interjacent[-1] + 1].location.start
            donor_end = anno.coordinate_genomic_to_gene(genomic_donor_end - 1, self.gene_id) + 1
        else:
            genomic_end = tx_model.exon[interjacent[-1]].location.end
            start = anno.coordinate_genomic_to_gene(genomic_end - 1, self.gene_id)
            genomic_start = tx_model.exon[interjacent[0]].location.start
            end = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id) + 1

            genomic_donor_end = self.long_exon_end
            donor_start = anno.coordinate_genomic_to_gene(genomic_donor_end - 1, self.gene_id)

            to_insert_full_exon = interjacent[0] == 0 \
                or tx_model.exon[interjacent[0] - 1].location.end <= self.long_exon_start

            if to_insert_full_exon:
                genomic_donor_start = self.long_exon_start
            else:
                genomic_donor_start = tx_model.exon[interjacent[0] - 1].location.end
            donor_end = anno.coordinate_genomic_to_gene(genomic_end, self.gene_id) + 1

        ref = str(gene_seq.seq[start])
        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

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

    def create_long_insertion(self, alignment:A3SSTranscriptAlignment,
            anno:gtf.GenomicAnnotation, gene_seq:dna.DNASeqRecord, var_id:str
            ) -> seqvar.VariantRecord:
        """ """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = alignment.tx_model

        if strand == 1:
            insert_pos_genomic = self.flanking_exon_end - 1
            insert_pos = anno.coordinate_genomic_to_gene(insert_pos_genomic, self.gene_id)
            donor_start = self.long_exon_start
            donor_start = anno.coordinate_genomic_to_gene(donor_start, self.gene_id)
            if alignment.upstream_index < len(tx_model.exon) - 1 \
                    and tx_model.exon[alignment.upstream_index + 1].location.start \
                        < self.long_exon_end:
                donor_end = tx_model.exon[alignment.upstream_index + 1].location.start
            else:
                donor_end = self.long_exon_end
            donor_end = anno.coordinate_genomic_to_gene(donor_end - 1, self.gene_id) + 1
        else:
            insert_pos_genomic = self.flanking_exon_start
            insert_pos = anno.coordinate_genomic_to_gene(insert_pos_genomic, self.gene_id)
            donor_start = self.long_exon_end
            donor_start = anno.coordinate_genomic_to_gene(donor_start - 1, self.gene_id)
            if alignment.upstream_index > 0 \
                    and tx_model.exon[alignment.upstream_index - 1].location.end \
                        > self.long_exon_start:
                donor_end = tx_model.exon[alignment.upstream_index - 1].location.end
            else:
                donor_end = self.long_exon_start
            donor_end = anno.coordinate_genomic_to_gene(donor_end, self.gene_id) + 1
        ref = str(gene_seq.seq[insert_pos])
        genomic_position = f'{chrom}:{insert_pos_genomic + 1}:{insert_pos_genomic + 2}'

        location = FeatureLocation(
            seqname=self.gene_id,
            start=insert_pos,
            end=insert_pos + 1
        )
        alt = '<Ins>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': donor_start,
            'END': donor_end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Insertion'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def create_short_deletion(self, alignment:A3SSTranscriptAlignment,
            spanning:int, interjacent:List[int], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = alignment.tx_model

        if strand == 1:
            if interjacent:
                genomic_start = tx_model.exon[interjacent[0]].location.start
            else:
                genomic_start = tx_model.exon[spanning].location.start
            start = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id)
            genomic_end = min(self.short_exon_start, tx_model.exon[spanning].location.end)
            end = anno.coordinate_gene_to_genomic(end - 1, self.gene_id) + 1
        else:
            if interjacent:
                genomic_end = tx_model.exon[interjacent[-1]].location.end
            else:
                genomic_end = tx_model.exon[spanning].location.end
            start = anno.coordinate_gene_to_genomic(genomic_end - 1, self.gene_id)
            genomic_start = max(self.short_exon_end, tx_model.exon[spanning].location.start)
            end = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id) + 1
        ref = str(gene_seq.seq[start])

        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

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

    def create_short_substitution(self, alignment:A3SSTranscriptAlignment,
            interjacent:List[int], anno:gtf.GenomicAnnotation,
            gene_seq:dna.DNASeqRecord, var_id:str) -> seqvar.VariantRecord:
        """ """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = alignment.tx_model

        if strand == 1:
            genomic_start = tx_model.exon[interjacent[0]].location.start
            start = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id)
            genomic_end = tx_model.exon[interjacent[-1]].location.end
            end = anno.coordinate_genomic_to_gene(genomic_end - 1, self.gene_id) + 1

            genomic_donor_start = self.short_exon_start
            donor_start = anno.coordinate_genomic_to_gene(genomic_donor_start, self.gene_id)
            # if the transcript does not have a downstream exon, the entire
            # short exon will be inserted.
            to_insert_full_exon = interjacent[-1] + 1 == len(tx_model.exon) \
                or tx_model.exon[interjacent[-1] + 1].location.start >= self.short_exon_end

            if to_insert_full_exon:
                genomic_donor_end = self.short_exon_end
            else:
                genomic_donor_end = tx_model.exon[interjacent[-1] + 1].location.start
            donor_end = anno.coordinate_genomic_to_gene(genomic_donor_end - 1, self.gene_id) + 1
        else:
            genomic_end = tx_model.exon[interjacent[-1]].location.end
            start = anno.coordinate_genomic_to_gene(genomic_end - 1, self.gene_id)
            genomic_start = tx_model.exon[interjacent[0]].location.start
            end = anno.coordinate_genomic_to_gene(genomic_start, self.gene_id) + 1

            genomic_donor_end = self.short_exon_end
            donor_start = anno.coordinate_genomic_to_gene(genomic_donor_end - 1, self.gene_id)

            to_insert_full_exon = interjacent[0] == 0 \
                or tx_model.exon[interjacent[0] - 1].location.end <= self.short_exon_start

            if to_insert_full_exon:
                genomic_donor_start = self.short_exon_start
            else:
                genomic_donor_start = tx_model.exon[interjacent[0] - 1].location.end
            donor_end = anno.coordinate_genomic_to_gene(genomic_end, self.gene_id) + 1

        ref = str(gene_seq.seq[start])
        genomic_position = f'{chrom}:{genomic_start + 1}:{genomic_end}'

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

    def create_short_insertion(self, alignment:A3SSTranscriptAlignment,
            anno:gtf.GenomicAnnotation, gene_seq:dna.DNASeqRecord, var_id:str
            ) -> seqvar.VariantRecord:
        """ """
        gene_model = anno.genes[self.gene_id]
        tx_id = alignment.tx_model.transcript_id
        chrom = gene_model.chrom
        strand = gene_model.strand
        tx_model = alignment.tx_model

        if strand == 1:
            insert_pos_genomic = self.flanking_exon_end - 1
            insert_pos = anno.coordinate_genomic_to_gene(insert_pos_genomic, self.gene_id)
            donor_start = self.short_exon_start
            donor_start = anno.coordinate_genomic_to_gene(donor_start, self.gene_id)
            if alignment.upstream_index < len(tx_model.exon) - 1 \
                    and tx_model.exon[alignment.upstream_index + 1].location.start \
                        < self.short_exon_end:
                donor_end = tx_model.exon[alignment.upstream_index + 1].location.start
            else:
                donor_end = self.short_exon_end
            donor_end = anno.coordinate_genomic_to_gene(donor_end - 1, self.gene_id) + 1
        else:
            insert_pos_genomic = self.flanking_exon_start
            insert_pos = anno.coordinate_genomic_to_gene(insert_pos_genomic, self.gene_id)
            donor_start = self.short_exon_end
            donor_start = anno.coordinate_genomic_to_gene(donor_start - 1, self.gene_id)
            if alignment.upstream_index > 0 \
                    and tx_model.exon[alignment.upstream_index - 1].location.end \
                        > self.short_exon_start:
                donor_end = tx_model.exon[alignment.upstream_index - 1].location.end
            else:
                donor_end = self.short_exon_start
            donor_end = anno.coordinate_genomic_to_gene(donor_end, self.gene_id) + 1
        ref = str(gene_seq.seq[insert_pos])
        genomic_position = f'{chrom}:{insert_pos_genomic + 1}:{insert_pos_genomic + 2}'

        location = FeatureLocation(
            seqname=self.gene_id,
            start=insert_pos,
            end=insert_pos + 1
        )
        alt = '<Ins>'
        attrs = {
            'TRANSCRIPT_ID': tx_id,
            'START': donor_start,
            'END': donor_end,
            'GENE_SYMBOL': gene_model.gene_name,
            'GENOMIC_POSITION': genomic_position
        }
        _type = 'Insertion'
        return seqvar.VariantRecord(location, ref, alt, _type, var_id, attrs)

    def convert_to_variant_records(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict, min_ijc:int, min_sjc:int
            ) -> List[seqvar.VariantRecord]:
        variants = []
        if not self.has_novel_splicing(anno):
            return variants

        gene_model = anno.genes[self.gene_id]
        tx_ids = gene_model.transcripts
        chrom = gene_model.location.seqname
        gene_seq = gene_model.get_gene_sequence(genome[chrom])

        tx_alignments:Dict[str, A3SSTranscriptAlignment] = {}

        for tx_id in tx_ids:
            tx_model = anno.transcripts[tx_id]

            alignment = self.align_to_transcript(tx_model)
            if alignment:
                tx_alignments[tx_id] = alignment

        var_id = self.create_variant_id(anno)

        for tx_id, alignment in tx_alignments.items():
            short_interjacent = alignment.get_short_interjacent_exons()
            long_interjacent = alignment.get_long_interjacent_exons()
            short_spanning = alignment.get_short_spanning_exon()
            long_spanning = alignment.get_long_spanning_exon()

            if alignment.long_index == -1:
                if long_spanning > -1:
                    record = self.create_long_deletion(
                        alignment, long_spanning, long_interjacent,
                        anno, gene_seq, var_id
                    )
                elif len(long_interjacent) > 0:
                    record = self.create_long_substitution(
                        alignment, long_interjacent,
                        anno, gene_seq, var_id
                    )
                else:
                    record = self.create_long_insertion(
                        alignment, anno, gene_seq, var_id
                    )
                variants.append(record)

            if alignment.short_index == -1:
                if short_spanning > -1:
                    record = self.create_short_deletion(
                        alignment, short_spanning, short_interjacent,
                        anno, gene_seq, var_id
                    )
                elif len(short_interjacent) > 0:
                    record = self.create_short_substitution(
                        alignment, short_interjacent,
                        anno, gene_seq, var_id
                    )
                else:
                    record = self.create_short_insertion(
                        alignment, anno, gene_seq, var_id
                    )
                variants.append(record)
        return variants
