""" Test module for VariantRecord """
import unittest
from test.unit import create_variant, create_genomic_annotation, \
    create_dna_record_dict, create_variants
from moPepGen import seqvar
from moPepGen.dna.DNASeqRecord import DNASeqRecord


GENOME_DATA = {
    'chr1': 'ATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCT'
}
ANNOTATION_ATTRS = [
    {
        'gene_id': 'ENSG0001'
    },{
        'transcript_id': 'ENST0001',
        'gene_id': 'ENSG0001',
        'protein_id': 'ENSP0001'
    }
]
ANNOTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 40, ANNOTATION_ATTRS[0]),
        'transcripts': ['ENST0001']
    }],
    'transcripts': [{
        'transcript_id': ANNOTATION_ATTRS[1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (5, 35, ANNOTATION_ATTRS[1]),
        'exon': [
            (5, 12, ANNOTATION_ATTRS[1]),
            (17, 23, ANNOTATION_ATTRS[1]),
            (27, 35, ANNOTATION_ATTRS[1])
        ]
    }]
}

class TestVariantRecord(unittest.TestCase):
    """ Test cases for VariantRecord """
    def test_to_transcript_variant_main(self):
        """ Test variant that contains an intron """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        variant_data = [
            10, 11, 'C', 'G', 'SNV', 'SNV-10-C-G',
            {'TRANSCRIPT_ID': 'ENST0001'}, 'ENSG0001'
        ]
        variant = create_variant(*variant_data)
        variant_tx = variant.to_transcript_variant(anno, genome, 'ENST0001')
        self.assertEqual(variant_tx.location.start, 5)
        self.assertEqual(variant_tx.location.end, 6)
        self.assertEqual(variant_tx.ref, 'C')
        self.assertEqual(variant_tx.alt, 'G')
        self.assertEqual(variant_tx.location.seqname, 'ENST0001')

    def test_to_transcript_variant_with_intron(self):
        """ Test variant that contains an intron """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        variant_data = [
            8, 20, 'GTCCTTCTGCCT', 'G', 'INDEL', 'INDEL-8-XXX',
            {'TRANSCRIPT_ID': 'ENST0001'}, 'ENSG0001'
        ]
        variant = create_variant(*variant_data)
        variant_tx = variant.to_transcript_variant(anno, genome, 'ENST0001')
        self.assertEqual(variant_tx.ref, 'GTCCCCT')

    def test_is_spanning_over_splicing_site(self):
        """ Check if the variant is spanning over a splicing site """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        variant_data = [
            21, 24, 'CCC', 'C', 'INDEL', 'INDEL-20-CCC-C',
            {'TRANSCRIPT_ID': 'ENST0001'}, 'ENSG0001'
        ]
        variant = create_variant(*variant_data)
        self.assertTrue(variant.is_spanning_over_splicing_site(anno, 'ENST0001'))

    def test_to_end_inclusion(self):
        """ Test that variant record is converted to end inclusion """
        tx_id = 'ENST0001'
        variant_data = [
            5, 6, 'C', 'CAAAA', 'INDEL', 'INDEL-5-C-CAAAA',
            {'TRANSCRIPT_ID': tx_id}, 'ENSG0001'
        ]
        variant = create_variant(*variant_data)
        seq = DNASeqRecord('C' * 10, id=tx_id, name=tx_id, description=tx_id)
        variant.to_end_inclusion(seq)
        self.assertEqual(variant.ref, 'C')
        self.assertEqual(variant.alt, 'AAAAC')
        self.assertEqual(variant.location.start, 6)

        # Insertion
        variant2 = create_variant(*variant_data)
        variant2.alt = '<Insertion>'
        variant2.type = 'Insertion'
        seq = DNASeqRecord('C' * 10, id=tx_id, name=tx_id, description=tx_id)
        variant2.to_end_inclusion(seq)
        self.assertEqual(variant2.ref, 'C')
        self.assertEqual(variant2.location.start, 6)
        self.assertTrue(variant2.is_end_inclusion())

    def test_to_end_inclusion_fail(self):
        """ Should fail when converting non insertion variants to end inclusion """
        tx_id = 'ENST0001'
        variant_data = [
            5, 6, 'C', '<Deletion>', 'Deletion', 'RI-5',
            {'TRANSCRIPT_ID': tx_id}, 'ENSG0001'
        ]
        variant = create_variant(*variant_data)
        seq = DNASeqRecord('C' * 10, id=tx_id, name=tx_id, description=tx_id)
        with self.assertRaises(ValueError):
            variant.to_end_inclusion(seq)

class TestTranscriptionalVariantSeries(unittest.TestCase):
    """ Test cases for VariantRecordSeries """
    def test_highest_hypermutated_region_complexity(self):
        """ Test calculating the highest hypermutated region complexity """
        var_data = {
            (3, 4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (4, 5, 'G', 'A', 'SNV', '', None, 'ENST0001.1'),
            (5, 6, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        variants = create_variants(var_data)
        series = seqvar.TranscriptionalVariantSeries(transcriptional=variants)
        series.sort()
        x = series.get_highest_hypermutated_region_complexity()
        self.assertEqual(x, 3)

        var_data = {
            (3, 4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (4, 5, 'G', 'A', 'SNV', '', None, 'ENST0001.1'),
            (5, 6, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (20, 21, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        variants = create_variants(var_data)
        series = seqvar.TranscriptionalVariantSeries(transcriptional=variants)
        series.sort()
        x = series.get_highest_hypermutated_region_complexity()
        self.assertEqual(x, 3)
