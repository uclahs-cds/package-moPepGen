""" Test module for reditools parser """
import copy
import unittest
from test.unit.test_vep_parser import ANNOTATION_ATTRS
from test.unit import create_genomic_annotation
from moPepGen import parser


ANNOTATION_ATTRS = [
    {
        'gene_id': 'ENSG0001',
        'gene_name': 'SYMBO'
    }, {
        'transcript_id': 'ENST0001',
        'gene_id': 'ENSG0001',
        'protein_id': 'ENSP0001',
        'gene_name': 'SYMBO'
    }
]
ANNOTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 400, ANNOTATION_ATTRS[0]),
        'transcripts': ['ENST0001']
    }],
    'transcripts': [{
        'transcript_id': ANNOTATION_ATTRS[1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (150, 350, ANNOTATION_ATTRS[1]),
        'exon': [
            (150, 200, ANNOTATION_ATTRS[1]),
            (225, 275, ANNOTATION_ATTRS[1]),
            (300, 350, ANNOTATION_ATTRS[1])
        ]
    }]
}

class TestREDItoolsParser(unittest.TestCase):
    """ Test case for REDItoolsParser """
    def test_reditools_record_to_variant_case1(self):
        """ strand positive """
        anno = create_genomic_annotation(ANNOTATION_DATA)

        # first exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=175,
            reference='C',
            strand=1,
            coverage_q=15,
            mean_quality=40.47,
            base_count=[0, 5, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            g_coverage_q=20,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(
            anno=anno, min_coverage_alt=3, min_frequency_alt=0.1,
            min_coverage_dna=10
        )
        self.assertEqual(variants[0].location.start, 174)
        self.assertEqual(variants[0].location.end, 175)
        self.assertEqual(variants[0].ref, 'C')
        self.assertEqual(variants[0].alt, 'T')

        # second exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=250,
            reference='C',
            strand=1,
            coverage_q=15,
            mean_quality=40.47,
            base_count=[0, 5, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            g_coverage_q=20,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(
            anno=anno, min_coverage_alt=3, min_frequency_alt=0.1,
            min_coverage_dna=10
        )
        self.assertEqual(variants[0].location.start, 249)
        self.assertEqual(variants[0].location.end, 250)
        self.assertEqual(variants[0].ref, 'C')
        self.assertEqual(variants[0].alt, 'T')

    def test_reditools_record_to_variant_case2(self):
        """ strand negative """
        annotation_data = copy.deepcopy(ANNOTATION_DATA)
        annotation_data['genes'][0]['strand'] = -1
        annotation_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(annotation_data)

        # first exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=325,
            reference='C',
            strand=-1,
            coverage_q=15,
            mean_quality=40.47,
            base_count=[0, 5, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            g_coverage_q=20,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(
            anno=anno, min_coverage_alt=3, min_frequency_alt=0.1,
            min_coverage_dna=10
        )
        self.assertEqual(variants[0].location.start, 75)
        self.assertEqual(variants[0].location.end, 76)
        self.assertEqual(variants[0].ref, 'G')
        self.assertEqual(variants[0].alt, 'A')

        # second exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=250,
            reference='C',
            strand=-1,
            coverage_q=15,
            mean_quality=40.47,
            base_count=[0, 5, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            g_coverage_q=20,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(
            anno=anno, min_coverage_alt=3, min_frequency_alt=0.1,
            min_coverage_dna=10
        )
        self.assertEqual(variants[0].location.start, 150)
        self.assertEqual(variants[0].location.end, 151)
        self.assertEqual(variants[0].ref, 'G')
        self.assertEqual(variants[0].alt, 'A')
