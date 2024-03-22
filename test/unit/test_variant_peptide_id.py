""" Test module for VariantPeptideIdentifier """
from typing import List
import unittest
from test.unit import create_variants
from moPepGen import aa
from moPepGen.aa import VariantPeptideIdentifier as pi


class TestVaraintPeptideIdentifier(unittest.TestCase):
    """ Test cases for VariantPeptideidentifier """
    def test_create_variant_id_single_snv(self):
        """ single SNV variant """
        variant_data = [
            (101, 102, 'A', 'T', 'SNV', 'SNV-101-A-T', None, 'ENST0001')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        self.assertEqual(str(peptide_id), 'ENST0001|SNV-101-A-T|1')

    def test_create_variant_id_multiple_snv(self):
        """ multiple snv """
        variant_data = [
            (101, 102, 'A', 'T', 'SNV', 'SNV-101-A-T', None, 'ENST0001'),
            (111, 112, 'A', 'T', 'SNV', 'SNV-111-A-T', None, 'ENST0001')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        self.assertEqual(str(peptide_id), 'ENST0001|SNV-101-A-T|SNV-111-A-T|1')

    def test_create_variant_id_single_fusion(self):
        """ single fusion variant """
        variant_data = [
            (101, 102, 'A', 'T', 'Fusion', 'FUSION-ENSG0001:101-ENSG0002:201', {
                'GENE_ID': 'ENSG0001',
                'ACCEPTER_GENE_ID': 'ENSG0002',
                'ACCEPTER_TRANSCRIPT_ID': 'ENST0002'
            }, 'ENST0001')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        self.assertEqual(str(peptide_id), 'FUSION-ENSG0001:101-ENSG0002:201|1')

    def test_create_variant_id_fusion_and_first_snv(self):
        """ fusion & 1-snv """
        variant_data = [
            (101, 102, 'A', 'T', 'Fusion', 'FUSION-ENSG0001:101-ENSG0002:201', {
                'GENE_ID': 'ENSG0001',
                'ACCEPTER_GENE_ID': 'ENSG0002',
                'ACCEPTER_TRANSCRIPT_ID': 'ENST0002'
            }, 'ENST0001'),
            (98, 99, 'A', 'T', 'SNV', 'SNV-98-A-T', None, 'ENST0001')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        expected = 'FUSION-ENSG0001:101-ENSG0002:201|1-SNV-98-A-T|1'
        self.assertEqual(str(peptide_id), expected)

    def test_create_variant_id_fusion_and_second_snv(self):
        """ fusion & 2-snv """
        variant_data = [
            (101, 102, 'A', 'T', 'Fusion', 'FUSION-ENSG0001:101-ENSG0002:201', {
                'GENE_ID': 'ENSG0001',
                'ACCEPTER_GENE_ID': 'ENSG0002',
                'ACCEPTER_TRANSCRIPT_ID': 'ENST0002'
            }, 'ENST0001'),
            (210, 211, 'A', 'T', 'SNV', 'SNV-210-A-T', None, 'ENST0002')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        expected = 'FUSION-ENSG0001:101-ENSG0002:201|2-SNV-210-A-T|1'
        self.assertEqual(str(peptide_id), expected)

    def test_create_variant_id_fusion_and_first_and_second_snv(self):
        """ fusion & 1-snv & 2-snv """
        variant_data = [
            (101, 102, 'A', 'T', 'Fusion', 'FUSION-ENSG0001:101-ENSG0002:201', {
                'GENE_ID': 'ENSG0001',
                'ACCEPTER_GENE_ID': 'ENSG0002',
                'ACCEPTER_TRANSCRIPT_ID': 'ENST0002'
            }, 'ENST0001'),
            (98, 99, 'A', 'T', 'SNV', 'SNV-98-A-T', None, 'ENST0001'),
            (210, 211, 'A', 'T', 'SNV', 'SNV-210-A-T', None, 'ENSG0002'),
            (215, 216, 'A', 'T', 'SNV', 'SNV-215-A-T', None, 'ENST0002')
        ]
        variants = create_variants(variant_data)
        variants[2].attrs['TRANSCRIPT_ID'] = 'ENST0002'
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        expected = 'FUSION-ENSG0001:101-ENSG0002:201|1-SNV-98-A-T|' +\
            '2-SNV-210-A-T|2-SNV-215-A-T|1'
        self.assertEqual(str(peptide_id), expected)

    def test_create_variant_id_single_circ_rna(self):
        """ single circRNA variant """
        variant_data = [
            (101, 102, 'A', '<circRNA>', 'circRNA', 'CIRC-ENSG0001-E2-E3', None, 'ENSG0001')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=1)
        self.assertEqual(str(peptide_id), 'CIRC-ENSG0001-E2-E3|1')

    def test_create_variant_id_single_circ_rna_and_snv(self):
        """ circRNA & snv """
        variant_data = [
            (101, 102, 'A', '<circRNA>', 'circRNA', 'CIRC-ENSG0001-E2-E3', None, 'ENSG0001'),
            (111, 112, 'A', 'T', 'SNV', 'SNV-111-A-T', None, 'ENST0001')
        ]
        variants = create_variants(variant_data)
        peptide_id = aa.create_variant_peptide_id('ENST0001', variants, index=2)
        self.assertEqual(str(peptide_id), 'CIRC-ENSG0001-E2-E3|SNV-111-A-T|2')

    def test_create_variant_id_orf(self):
        """ orf """
        peptide_id = aa.create_variant_peptide_id('ENST0001', [], 'ORF1', 1)
        self.assertEqual(str(peptide_id), 'ENST0001|ORF1|1')

    def test_parse_variant_id_single_snv(self):
        """ parse variant from single snv """
        coding_txs = {}
        label = 'ENST0001|SNV-101-A-T|1'
        peptide_ids:List[pi.BaseVariantPeptideIdentifier]
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.BaseVariantPeptideIdentifier)
        self.assertEqual(peptide_ids[0].transcript_id, 'ENST0001')
        self.assertEqual(set(peptide_ids[0].variant_ids), {'SNV-101-A-T'})
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_multiple_snv(self):
        """ parse variant from multiple snv """
        coding_txs = {}
        label = 'ENST0001|SNV-101-A-T|SNV-111-A-T|1'
        peptide_ids:List[pi.BaseVariantPeptideIdentifier]
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.BaseVariantPeptideIdentifier)
        self.assertEqual(peptide_ids[0].transcript_id, 'ENST0001')
        self.assertEqual(set(peptide_ids[0].variant_ids), {'SNV-101-A-T', 'SNV-111-A-T'})
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_single_fusion(self):
        """ parse variant from single fusion """
        coding_txs = {}
        label = 'FUSION-ENSG0001:101-ENSG0002:201|1'
        peptide_ids:List[pi.FusionVariantPeptideIdentifier]
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.FusionVariantPeptideIdentifier)
        self.assertEqual(peptide_ids[0].fusion_id, 'FUSION-ENSG0001:101-ENSG0002:201')
        self.assertEqual(len(peptide_ids[0].first_variants), 0)
        self.assertEqual(len(peptide_ids[0].second_variants), 0)
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_fusion_and_first_variant(self):
        """ parse variant from single fusion """
        coding_txs = {}
        label = 'FUSION-ENSG0001:101-ENSG0002:201|1-SNV-98-A-T|1'
        peptide_ids:List[pi.FusionVariantPeptideIdentifier]
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.FusionVariantPeptideIdentifier)
        self.assertEqual(peptide_ids[0].fusion_id, 'FUSION-ENSG0001:101-ENSG0002:201')
        self.assertEqual(set(peptide_ids[0].first_variants), {'SNV-98-A-T'})
        self.assertEqual(len(peptide_ids[0].second_variants), 0)
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_fusion_and_second_variant(self):
        """ parse variant from single fusion """
        coding_txs = {}
        label = 'FUSION-ENSG0001:101-ENSG0002:201|2-SNV-210-A-T|1'
        peptide_ids:List[pi.FusionVariantPeptideIdentifier]
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.FusionVariantPeptideIdentifier)
        self.assertEqual(peptide_ids[0].fusion_id, 'FUSION-ENSG0001:101-ENSG0002:201')
        self.assertEqual(len(peptide_ids[0].first_variants), 0)
        self.assertEqual(set(peptide_ids[0].second_variants), {'SNV-210-A-T'})
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_fusion_and_first_and_second_variant(self):
        """ parse variant from single fusion """
        coding_txs = {}
        label = 'FUSION-ENSG0001:101-ENSG0002:201|1-SNV-98-A-T|' +\
            '2-SNV-210-A-T|2-SNV-215-A-T|1'
        peptide_ids:List[pi.FusionVariantPeptideIdentifier]
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.FusionVariantPeptideIdentifier)
        self.assertEqual(peptide_ids[0].fusion_id, 'FUSION-ENSG0001:101-ENSG0002:201')
        self.assertEqual(set(peptide_ids[0].first_variants), {'SNV-98-A-T'})
        self.assertEqual(set(peptide_ids[0].second_variants), {'SNV-210-A-T','SNV-215-A-T'})
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_single_circ_rna(self):
        """ parse variant from single fusion """
        coding_txs = {}
        label = 'CIRC-ENSG0001-E2-E3|1'
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.CircRNAVariantPeptideIdentifier)
        peptide_ids:List[pi.CircRNAVariantPeptideIdentifier]
        self.assertEqual(peptide_ids[0].circ_rna_id, 'CIRC-ENSG0001-E2-E3')
        self.assertEqual(len(peptide_ids[0].variant_ids), 0)
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_circ_rna_and_snv(self):
        """ parse variant from single fusion """
        coding_txs = {}
        label = 'CIRC-ENSG0001-E2-E3|SNV-111-A-T|2'
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.CircRNAVariantPeptideIdentifier)
        peptide_ids:List[pi.CircRNAVariantPeptideIdentifier]
        self.assertEqual(peptide_ids[0].circ_rna_id, 'CIRC-ENSG0001-E2-E3')
        self.assertEqual(set(peptide_ids[0].variant_ids), {'SNV-111-A-T'})
        self.assertEqual(str(peptide_ids[0]), label)

    def test_parse_variant_id_orf(self):
        """ parse variant with orf """
        coding_txs = {}
        label = 'ENST0001|ENSG0001|ORF1|1'
        peptide_ids = aa.parse_variant_peptide_id(label, coding_txs)
        self.assertEqual(len(peptide_ids), 1)
        self.assertIsInstance(peptide_ids[0], pi.NovelORFPeptideIdentifier)
        peptide_ids:List[pi.NovelORFPeptideIdentifier]
        self.assertEqual(peptide_ids[0].transcript_id, 'ENST0001')
        self.assertEqual(peptide_ids[0].orf_id, 'ORF1')
        self.assertEqual(str(peptide_ids[0]), label)
