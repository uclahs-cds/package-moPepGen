""" Test module for ArribaParser """
import copy
import unittest
from test.unit import create_genomic_annotation, create_dna_record_dict, \
    GENOME_DATA, ANNOTATION_DATA
from moPepGen.parser.ArribaParser import ArribaConfidence, ArribaRecord


class TestArribaParser(unittest.TestCase):
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
        self.assertEqual(variants[0].location.start, 23)
        tx_variant = variants[0].to_transcript_variant(anno, genome)
        self.assertEqual(tx_variant.location.start, 13)
