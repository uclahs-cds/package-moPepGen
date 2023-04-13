""" Test case for RMATRecord """
import copy
import unittest
from test.unit import create_genomic_annotation, create_dna_record_dict
from moPepGen.parser import RMATSParser


GENOME_DATA = {
    'chr1': 'ATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCT'
}
ANNOTATION_ATTRS = [
    {
        'gene_id': 'ENSG0001',
        'gene_name': 'SYMBO'
    },{
        'transcript_id': 'ENST0001.1',
        'gene_id': 'ENSG0001',
        'protein_id': 'ENSP0001',
        'gene_name': 'SYMBO'
    }
]
ANNOTATION_DATA  = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 50, ANNOTATION_ATTRS[0]),
        'transcripts': ['ENST0001.1']
    }],
    'transcripts': [{
        # seq: CTGGT CCCCT ATGGG TCCTT C
        'transcript_id': ANNOTATION_ATTRS[1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (5, 48, ANNOTATION_ATTRS[1]),
        'exon': [
            ( 5, 12, ANNOTATION_ATTRS[1]),
            (17, 23, ANNOTATION_ATTRS[1]),
            (27, 35, ANNOTATION_ATTRS[1]),
            (40, 48, ANNOTATION_ATTRS[1])
        ]
    }]
}

def create_a3ss_pos() -> RMATSParser.A5SSRecord:
    """ A3SS with positive strand """
    return RMATSParser.A3SSRecord(
        gene_id='ENSG0001',
        gene_symbol='CRAP',
        chrom='chr1',
        long_exon_start=17,
        long_exon_end=23,
        short_exon_start=20,
        short_exon_end=23,
        flanking_exon_start=5,
        flanking_exon_end=12,
        ijc_sample_1=1,
        sjc_sample_1=1,
        ijc_sample_2='',
        sjc_sample_2='',
        inc_form_len=300,
        skip_form_len=74,
        pvalue='NA',
        fdr='NA'
    )


def create_a3ss_neg() -> RMATSParser.A5SSRecord:
    """ Create A3SS with negative strand """
    return RMATSParser.A3SSRecord(
        gene_id='ENSG0001',
        gene_symbol='CRAP',
        chrom='chr1',
        long_exon_start=17,
        long_exon_end=23,
        short_exon_start=17,
        short_exon_end=20,
        flanking_exon_start=27,
        flanking_exon_end=35,
        ijc_sample_1=1,
        sjc_sample_1=1,
        ijc_sample_2='',
        sjc_sample_2='',
        inc_form_len=300,
        skip_form_len=74,
        pvalue='NA',
        fdr='NA'
    )

class TestRMATSRecord(unittest.TestCase):
    """ Tset cases for RMATSRecord """
    def test_se_record_pos_strand(self):
        """ Test SERecord """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.SERecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            exon_start=17,
            exon_end=23,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=27,
            downstream_exon_end=35,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        gene_seq = anno.genes[gene_id].get_gene_sequence(genome[chrom])
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 6)
        self.assertEqual(var_records[0].ref, str(gene_seq.seq[17]))
        self.assertEqual(var_records[0].attrs['START'], 17)
        self.assertEqual(var_records[0].attrs['END'], 23)

    def test_se_record_neg_strand(self):
        """ Test SERecord on - strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.SERecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            exon_start=17,
            exon_end=23,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=27,
            downstream_exon_end=35,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        gene_seq = anno.genes[gene_id].get_gene_sequence(genome[chrom])
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 6)
        self.assertEqual(var_records[0].ref, str(gene_seq.seq[27]))
        self.assertEqual(var_records[0].attrs['START'], 27)
        self.assertEqual(var_records[0].attrs['END'], 33)

    def test_se_record_has_skipped_pos(self):
        """ Test SERecord """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        exons:list = anno_data['transcripts'][0]['exon']
        exons.pop(1)
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.SERecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            exon_start=17,
            exon_end=23,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=27,
            downstream_exon_end=35,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        gene_seq = anno.genes[gene_id].get_gene_sequence(genome[chrom])
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(var_records[0].location.start, 11)
        self.assertEqual(var_records[0].ref, str(gene_seq.seq[11]))
        self.assertEqual(var_records[0].type, 'Insertion')
        self.assertEqual(var_records[0].attrs['DONOR_START'], 17)
        self.assertEqual(var_records[0].attrs['DONOR_END'], 23)

    def test_se_record_has_skipped_neg(self):
        """ Test SERecord """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        exons:list = anno_data['transcripts'][0]['exon']
        exons.pop(1)
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.SERecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            exon_start=17,
            exon_end=23,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=27,
            downstream_exon_end=35,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        gene_seq = anno.genes[gene_id].get_gene_sequence(genome[chrom])
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(var_records[0].location.start, 22)
        self.assertEqual(var_records[0].ref, str(gene_seq.seq[22]))
        self.assertEqual(var_records[0].type, 'Insertion')
        self.assertEqual(var_records[0].attrs['DONOR_START'], 27)
        self.assertEqual(var_records[0].attrs['DONOR_END'], 33)

    def test_a5ss_record_pos_strand(self):
        """ Test A5SSRecord with pos strand """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.A5SSRecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            long_exon_start=17,
            long_exon_end=23,
            short_exon_start=17,
            short_exon_end=20,
            flanking_exon_start=27,
            flanking_exon_end=35,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 3)
        self.assertEqual(var_records[0].location.start, 20)
        self.assertEqual(var_records[0].location.end, 23)

    def test_a5ss_record_neg_strand(self):
        """ Test A5SSRecord with neg strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.A5SSRecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            long_exon_start=17,
            long_exon_end=23,
            short_exon_start=20,
            short_exon_end=23,
            flanking_exon_start=5,
            flanking_exon_end=12,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 3)
        self.assertEqual(var_records[0].location.start, 29)
        self.assertEqual(var_records[0].location.end, 32)

    def test_a5ss_record_has_short(self):
        """ Test that reference has the short version """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = RMATSParser.A5SSRecord(
            gene_id='ENSG0001',
            gene_symbol='CRAP',
            chrom='chr1',
            long_exon_start=17,
            long_exon_end=25,
            short_exon_start=17,
            short_exon_end=23,
            flanking_exon_start=27,
            flanking_exon_end=35,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(var_records[0].location.start, 22)
        self.assertEqual(var_records[0].attrs['DONOR_START'], 23)
        self.assertEqual(var_records[0].attrs['DONOR_END'], 25)

    def test_a3ss_pos_long(self):
        """ Test A5SSRecord with pos strand and long exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = create_a3ss_pos()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 3)
        self.assertEqual(var_records[0].location.start, 17)
        self.assertEqual(var_records[0].location.end, 20)

    def test_a3ss_pos_short(self):
        """ Test A5SSRecord with pos strand and short exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][1] = (20, 23, ANNOTATION_ATTRS[1])
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_pos()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(var_records[0].type, 'Insertion')
        self.assertEqual(var_records[0].attrs['START'], 17)
        self.assertEqual(var_records[0].attrs['END'], 20)

    def test_a3ss_pos_long_spanning(self):
        """ Test A5SSRecord with pos strand, long exon and spanning exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][1] = (15, 23, ANNOTATION_ATTRS[1])
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_pos()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertTrue(all(x.type == 'Deletion' for x in var_records))
        self.assertTrue(all(x.location.start == 15 for x in var_records))
        self.assertTrue({x.location.end for x in var_records} == {17, 20})

    def test_a3ss_pos_long_interjacent(self):
        """ Test A5SSRecord with pos strand, long exon and interjacent exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'].insert(1, (14, 15, ANNOTATION_ATTRS[1]))
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_pos()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertTrue(all(x.type == 'Deletion' for x in var_records))

    def test_a3ss_neg_long(self):
        """ Test A5SSRecord with neg strand, long exon """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_neg()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 3)
        self.assertEqual(var_records[0].location.start, 27)
        self.assertEqual(var_records[0].location.end, 30)

    def test_a3ss_neg_long_spanning(self):
        """ Test A5SSRecord with neg strand, long exon, and spanning exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][1] = (17, 25, ANNOTATION_ATTRS[1])
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)

        record = create_a3ss_neg()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertTrue(all(x.type == 'Deletion' for x in var_records))
        self.assertTrue(all(x.location.start == 25 for x in var_records))
        self.assertTrue({x.location.end for x in var_records} == {27, 30})

    def test_a3ss_neg_long_interjacent(self):
        """ Test A5SSRecord with neg strand, long exon, and an interjacent exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][2] = (30, 35, ANNOTATION_ATTRS[1])
        anno_data['transcripts'][0]['exon'].insert(2, (26, 27, ANNOTATION_ATTRS[1]))
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)

        record = create_a3ss_neg()
        record.flanking_exon_start = 30
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertTrue(all(x.type == 'Deletion' for x in var_records))

    def test_a3ss_pos_short_spanning(self):
        """ Test A5SSRecord with pos strand, short exon, and spanning exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][1] = (19, 23, ANNOTATION_ATTRS[1])
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_pos()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertEqual({x.type for x in var_records}, {'Insertion', 'Deletion'})

    def test_a3ss_pos_short_interjacent(self):
        """ Test A5SSRecord with positive strand, short exon and interjacent exons """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][1] = (20, 23, ANNOTATION_ATTRS[1])
        anno_data['transcripts'][0]['exon'].insert(1, (14, 15, ANNOTATION_ATTRS[1]))
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_pos()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertEqual({x.type for x in var_records}, {'Substitution', 'Deletion'})

    def test_a3ss_neg_short_spanning(self):
        """ Test A5SS with neg strand and spanning exon """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][1] = (17, 21, ANNOTATION_ATTRS[1])
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_neg()
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 2)
        self.assertEqual({x.type for x in var_records}, {'Insertion', 'Deletion'})

    def test_a3ss_neg_short_interjacent(self):
        """ Test A5SS with neg strand, short exon and interjacent exon """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['transcripts'][0]['exon'][2] = (30, 35, ANNOTATION_ATTRS[1])
        anno_data['transcripts'][0]['exon'].insert(2, (26, 27, ANNOTATION_ATTRS[1]))
        anno = create_genomic_annotation(anno_data)
        record = create_a3ss_pos()
        record.flanking_exon_start = 30
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertTrue(var_records[0].type, 'Substitution')

    def test_mxe_record_pos_strand(self):
        """ Test MXE with pos strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        exons:list = anno_data['transcripts'][0]['exon']
        exons.pop(1)
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.MXERecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            first_exon_start=17,
            first_exon_end=23,
            second_exon_start=27,
            second_exon_end=35,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=40,
            downstream_exon_end=48,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 8)
        self.assertEqual(var_records[0].location.start, 27)
        self.assertEqual(var_records[0].location.end, 35)

    def test_mxe_record_neg_strand(self):
        """ Test A5SSRecord with neg strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        exons:list = anno_data['transcripts'][0]['exon']
        exons.pop(1)
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.MXERecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            first_exon_start=17,
            first_exon_end=23,
            second_exon_start=27,
            second_exon_end=35,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=40,
            downstream_exon_end=48,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 8)
        self.assertEqual(var_records[0].location.start, 15)
        self.assertEqual(var_records[0].location.end, 23)

    def test_ri_record_pos_strand(self):
        """ Test RIRecord with pos strand """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.RIRecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            retained_intron_exon_start=12,
            retained_intron_exon_end=17,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=17,
            downstream_exon_end=23,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 1)
        self.assertEqual(var_records[0].attrs['DONOR_START'], 12)
        self.assertEqual(var_records[0].attrs['DONOR_END'], 17)

    def test_ri_record_neg_strand(self):
        """ Test A5SSRecord with neg strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.RIRecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            retained_intron_exon_start=12,
            retained_intron_exon_end=17,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=17,
            downstream_exon_end=23,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(len(var_records[0].location), 1)
        self.assertEqual(var_records[0].attrs['DONOR_START'], 33)
        self.assertEqual(var_records[0].attrs['DONOR_END'], 38)

    def test_ri_record_pos_spliced(self):
        """ Test RIRecord with pos strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        exons:list = anno_data['transcripts'][0]['exon']
        exon_1 = list(exons[0])
        exon_1[1] = 23
        exons[0] = tuple(exon_1)
        exons.pop(1)
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        gene_id = 'ENSG0001'
        chrom = 'chr1'
        record = RMATSParser.RIRecord(
            gene_id=gene_id,
            gene_symbol='CRAP',
            chrom=chrom,
            retained_intron_exon_start=12,
            retained_intron_exon_end=17,
            upstream_exon_start=5,
            upstream_exon_end=12,
            downstream_exon_start=17,
            downstream_exon_end=23,
            ijc_sample_1=1,
            sjc_sample_1=1,
            ijc_sample_2='',
            sjc_sample_2='',
            inc_form_len=300,
            skip_form_len=74,
            pvalue='NA',
            fdr='NA'
        )
        var_records = record.convert_to_variant_records(anno, genome, 1, 1)
        self.assertEqual(len(var_records), 1)
        self.assertEqual(var_records[0].location.start, 12)
        self.assertEqual(var_records[0].location.end, 17)
