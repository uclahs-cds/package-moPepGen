""" Test module for ArribaParser """
import copy
import unittest
from test.unit import create_genomic_annotation, create_dna_record_dict
from moPepGen.parser.ArribaParser import ArribaConfidence, ArribaRecord

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
        'transcript': (5, 35, ANNOTATION_ATTRS[1]['transcripts'][0]),
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
        fusion_record = ArribaRecord(
            gene1 = 'SYMBO1',
            gene2 = 'SYMBO2',
            strand1 = '+/+',
            strand2 = '+/+',
            breakpoint1 = 'chr1:21',
            breakpoint2 = 'chr1:81',
            site1 = '',
            site2 = '',
            _type = '',
            split_reads1 = '',
            split_reads2 = '',
            discordant_mates = 2,
            coverage1 = 2,
            coverage2 = 2,
            confidence = ArribaConfidence('high'),
            reading_frame = '',
            tags = '',
            retained_protein_domains = '',
            closest_genomic_breakpoint1 = '',
            closest_genomic_breakpoint2 = '',
            gene_id1 = 'ENSG0001',
            gene_id2 = 'ENSG0002',
            transcript_id1 = 'ENST0001.1',
            transcript_id2 = 'ENST0002.1',
            direction1 = 'downstream',
            direction2 = 'upstream',
            filters = [],
            fusion_transcript = '',
            peptide_sequence = '',
            read_identifiers = []
        )
        variants = fusion_record.convert_to_variant_records(anno, genome)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].ref, genome['chr1'][21])
        self.assertEqual(variants[0].attrs['ACCEPTER_POSITION'], 20)

    def test_convert_to_variant_negative_strand(self):
        """ Test convert to variant """
        genome = create_dna_record_dict(GENOME_DATA)
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)

        fusion_record = ArribaRecord(
            gene1 = 'SYMBO1',
            gene2 = 'SYMBO2',
            strand1 = '-/-',
            strand2 = '-/-',
            breakpoint1 = 'chr1:21',
            breakpoint2 = 'chr1:81',
            site1 = '',
            site2 = '',
            _type = '',
            split_reads1 = '',
            split_reads2 = '',
            discordant_mates = 2,
            coverage1 = 2,
            coverage2 = 2,
            confidence = ArribaConfidence('high'),
            reading_frame = '',
            tags = '',
            retained_protein_domains = '',
            closest_genomic_breakpoint1 = '',
            closest_genomic_breakpoint2 = '',
            gene_id1 = 'ENSG0001',
            gene_id2 = 'ENSG0002',
            transcript_id1 = 'ENST0001.1',
            transcript_id2 = 'ENST0002.1',
            direction1 = 'downstream',
            direction2 = 'upstream',
            filters = [],
            fusion_transcript = '',
            peptide_sequence = '',
            read_identifiers = []
        )
        variants = fusion_record.convert_to_variant_records(anno, genome)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].ref, genome['chr1'].seq[21:22].reverse_complement())
        self.assertEqual(variants[0].attrs['ACCEPTER_POSITION'], 20)

    def test_convert_to_variant_splice_site(self):
        """ Test convert to variant """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        fusion_record = ArribaRecord(
            gene1 = 'SYMBO1',
            gene2 = 'SYMBO2',
            strand1 = '+/+',
            strand2 = '+/+',
            breakpoint1 = 'chr1:23',
            breakpoint2 = 'chr1:78',
            site1 = '',
            site2 = '',
            _type = '',
            split_reads1 = '',
            split_reads2 = '',
            discordant_mates = 2,
            coverage1 = 2,
            coverage2 = 2,
            confidence = ArribaConfidence('high'),
            reading_frame = '',
            tags = '',
            retained_protein_domains = '',
            closest_genomic_breakpoint1 = '',
            closest_genomic_breakpoint2 = '',
            gene_id1 = 'ENSG0001',
            gene_id2 = 'ENSG0002',
            transcript_id1 = 'ENST0001.1',
            transcript_id2 = 'ENST0002.1',
            direction1 = 'downstream',
            direction2 = 'upstream',
            filters = [],
            fusion_transcript = '',
            peptide_sequence = '',
            read_identifiers = []
        )
        variants = fusion_record.convert_to_variant_records(anno, genome)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].location.start, 22)
        tx_variant = variants[0].to_transcript_variant(anno, genome)
        self.assertEqual(tx_variant.location.start, 13)
