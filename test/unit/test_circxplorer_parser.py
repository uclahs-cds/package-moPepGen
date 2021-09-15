""" Test module for CIRCexplorerParser """
import unittest
from test.unit import create_genomic_annotation
from test.unit.test_vep_parser import ANNOTATION_DATA
from moPepGen.parser.CIRCexplorerParser import CIRCexplorerKnownRecord
from moPepGen.circ import CircRNAModel


class TestCIRCexplorerParser(unittest.TestCase):
    """ Test cases fro CIRCexplorerParser """
    @staticmethod
    def create_base_record():
        """ create a base CIRCexplorerKnownRecord """
        return CIRCexplorerKnownRecord(
            chrom='chr1',
            start=None,
            end=None,
            name=None,
            score=None,
            strand='+',
            thick_start=None,
            thick_end=None,
            item_rgb=(0,0,0),
            exon_count=None,
            exon_sizes=None,
            exon_offsets=None,
            read_number=None,
            circ_type='circRNA',
            gene_name='ENSG0001',
            isoform_name='ENST0001.1',
            index=None,
            flank_intron=None
        )

    def test_to_convert_circ_rna(self):
        """ circRNA """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = self.create_base_record()
        record.start = 5
        record.end = 23
        record.exon_sizes = [7, 6]
        record.exon_offsets = [0, 12]
        circ_record = record.convert_to_circ_rna(anno)
        self.assertIsInstance(circ_record, CircRNAModel)
        self.assertEqual(circ_record.id, 'CIRC-ENSG0001-E1-E2')

    def test_to_convert_ci_rna(self):
        """ ciRNA """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = self.create_base_record()
        record.start = 12
        record.end = 17
        record.exon_sizes = [5]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        circ_record = record.convert_to_circ_rna(anno)
        self.assertIsInstance(circ_record, CircRNAModel)
        self.assertEqual(circ_record.id, 'CI-ENSG0001-I1')
