""" Test module for FusionCatcherParser """
import copy
import unittest
from test.unit import create_genomic_annotation, create_dna_record_dict
from Bio.Seq import Seq
from moPepGen.parser.FusionCatcherParser import FusionCatcherRecord

GENOME_DATA = {
    'chr1':
    'ATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCT'
    'CCTCCCAATAAAGTCGAATTTTGGAACCGAATTCCCTTTTTTCGGGAAAAGCTACTAGGG'
}
ANNOTATION_ATTRS = [
    {
        'gene': {
            'gene_id'  : 'ENSG0001',
            'gene_name': 'SYMBO1'
        },
        'transcripts': [{
            'transcript_id': 'ENST0001.1',
            'gene_id'      : 'ENSG0001',
            'protein_id'   : 'ENSP0001',
            'gene_name'    : 'SYMBO1'
        }]
    }, {
        'gene': {
            'gene_id'  : 'ENSG0002',
            'gene_name': 'SYMBO2'
        },
        'transcripts': [{
            'transcript_id': 'ENST0002.1',
            'gene_id'      : 'ENSG0002',
            'protein_id'   : 'ENSP0002',
            'gene_name'    : 'SYMBO2'
        }]
    }
]
ANNOTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene']['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 40, ANNOTATION_ATTRS[0]['gene']),
        'transcripts': ['ENST0001.1']
    }, {
        'gene_id': ANNOTATION_ATTRS[1]['gene']['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (60, 100, ANNOTATION_ATTRS[1]['gene']),
        'transcripts': ['ENST0002.1']
    }],
    'transcripts': [{
        'transcript_id': ANNOTATION_ATTRS[0]['transcripts'][0]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (5, 35, ANNOTATION_ATTRS[0]['transcripts'][0]),
        'exon': [
            (5,  12, ANNOTATION_ATTRS[0]['transcripts'][0]),
            (17, 23, ANNOTATION_ATTRS[0]['transcripts'][0]),
            (27, 35, ANNOTATION_ATTRS[0]['transcripts'][0])
        ]
    }, {
        'transcript_id': ANNOTATION_ATTRS[1]['transcripts'][0]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (65, 95, ANNOTATION_ATTRS[1]['transcripts'][0]),
        'exon': [
            (65, 72, ANNOTATION_ATTRS[1]['transcripts'][0]),
            (77, 83, ANNOTATION_ATTRS[1]['transcripts'][0]),
            (87, 95, ANNOTATION_ATTRS[1]['transcripts'][0])
        ]
    }]
}

class TestFusionCatcherParser(unittest.TestCase):
    """ Test cases for FusionCatcherParser """
    def test_convert_to_variant(self):
        """ Test convert to variant """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        fusion_record = FusionCatcherRecord(
            five_end_gene_symbol = 'SYMBO1',
            three_end_gene_symbol = 'SYMBOL2',
            fusion_descriptions = '',
            counts_of_common_mapping_reads =2,
            spanning_pairs = 2,
            spanning_unique_reads = 2,
            longest_anchor_found = 2,
            fusion_finding_method = '',
            five_end_breakpoint = '1:21:+',
            three_end_breakpoint = '1:81:+',
            five_end_gene_id = 'ENSG0001',
            three_end_gene_id = 'ENSG0002',
            five_end_exon_id = '',
            three_end_exon_id = '',
            fusion_sequence = 'CTTCTGCCTA*TTGGAACCGA',
            predicted_effect = ''
        )
        variants = fusion_record.convert_to_variant_records(anno, genome)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].ref, genome['chr1'][21])

    def test_convert_to_variant_negative_strand(self):
        """ Test convert to variant """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)

        fusion_record = FusionCatcherRecord(
            five_end_gene_symbol = 'SYMBO1',
            three_end_gene_symbol = 'SYMBOL2',
            fusion_descriptions = '',
            counts_of_common_mapping_reads =2,
            spanning_pairs = 2,
            spanning_unique_reads = 2,
            longest_anchor_found = 2,
            fusion_finding_method = '',
            five_end_breakpoint = '1:21:-',
            three_end_breakpoint = '1:81:+',
            five_end_gene_id = 'ENSG0001',
            three_end_gene_id = 'ENSG0002',
            five_end_exon_id = '',
            three_end_exon_id = '',
            fusion_sequence = 'CTTCTGCCTA*TTGGAACCGA',
            predicted_effect = ''
        )
        variants = fusion_record.convert_to_variant_records(anno, genome)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].ref, Seq(genome['chr1'][19]).reverse_complement())
