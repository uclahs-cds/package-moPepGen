""" Test module for VariantRecord """
import unittest
from test.unit import create_variant, create_genomic_annotation, \
    create_dna_record_dict


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
