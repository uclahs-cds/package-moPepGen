""" Test module for reditools parser """
import unittest
from test.unit import create_transcript_model
from moPepGen import parser


class TestREDItoolsParser(unittest.TestCase):
    """ Test case for REDItoolsParser """
    def test_reditools_record_to_variant_case1(self):
        """ strand positive """
        attributes = {
            'transcript_id': 'ENST0002',
            'gene_id': 'ENSG0002',
            'protein_id': 'ENSP0002'
        }
        data = {
            'chrom': 'chr22',
            'strand': 1,
            'transcript': (150, 350, attributes),
            'exon': [
                (150,200,attributes),
                (225, 275, attributes),
                (300, 350, attributes)
            ]
        }
        model = create_transcript_model(data)
        anno = {'ENST0001': model}

        # first exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=175,
            reference='C',
            strand=1,
            coverage_q30=15,
            mean_quality=40.47,
            base_count=[0, 2, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(anno)
        self.assertEqual(variants[0].location.start, 24)
        self.assertEqual(variants[0].location.end, 25)
        self.assertEqual(variants[0].ref, 'C')
        self.assertEqual(variants[0].alt, 'T')

        # second exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=250,
            reference='C',
            strand=1,
            coverage_q30=15,
            mean_quality=40.47,
            base_count=[0, 2, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(anno)
        self.assertEqual(variants[0].location.start, 74)
        self.assertEqual(variants[0].location.end, 75)
        self.assertEqual(variants[0].ref, 'C')
        self.assertEqual(variants[0].alt, 'T')

    def test_reditools_record_to_variant_case2(self):
        """ strand negative """
        attributes = {
            'transcript_id': 'ENST0002',
            'gene_id': 'ENSG0002',
            'protein_id': 'ENSP0002'
        }
        data = {
            'chrom': 'chr22',
            'strand': -1,
            'transcript': (150, 350, attributes),
            'exon': [
                (150,200,attributes),
                (225, 275, attributes),
                (300, 350, attributes)
            ]
        }
        model = create_transcript_model(data)
        anno = {'ENST0001': model}

        # first exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=325,
            reference='C',
            strand=-1,
            coverage_q30=15,
            mean_quality=40.47,
            base_count=[0, 2, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(anno)
        self.assertEqual(variants[0].location.start, 25)
        self.assertEqual(variants[0].location.end, 26)
        self.assertEqual(variants[0].ref, 'G')
        self.assertEqual(variants[0].alt, 'A')

        # second exon
        record = parser.REDItoolsRecord(
            region='chr22',
            position=250,
            reference='C',
            strand=-1,
            coverage_q30=15,
            mean_quality=40.47,
            base_count=[0, 2, 0, 13],
            all_subs=[('C', 'T')],
            frequency=0.87,
            transcript_id=[('ENST0001', 'transcript')]
        )
        variants = record.convert_to_variant_records(anno)
        self.assertEqual(variants[0].location.start, 75)
        self.assertEqual(variants[0].location.end, 76)
        self.assertEqual(variants[0].ref, 'G')
        self.assertEqual(variants[0].alt, 'A')
