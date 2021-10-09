""" Test module for CIRCexplorerParser """
import unittest
import copy
from test.unit import create_genomic_annotation
from test.unit.test_vep_parser import ANNOTATION_ATTRS, ANNOTATION_DATA
from moPepGen.parser.CIRCexplorerParser import CIRCexplorer2KnownRecord
from moPepGen.circ import CircRNAModel
from moPepGen import err



class TestCIRCexplorerParser(unittest.TestCase):
    """ Test cases fro CIRCexplorerParser """
    @staticmethod
    def create_base_record():
        """ create a base CIRCexplorerKnownRecord """
        return CIRCexplorer2KnownRecord(
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
        self.assertEqual(circ_record.id, 'CIRC-ENST0001.1-E1-E2')

    def test_to_convert_circ_rna_negative(self):
        """ circRNA on -strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)
        record = self.create_base_record()
        record.start = 5
        record.end = 23
        record.exon_sizes = [7, 6]
        record.exon_offsets = [0, 12]
        circ_record = record.convert_to_circ_rna(anno)
        self.assertIsInstance(circ_record, CircRNAModel)
        self.assertEqual(circ_record.id, 'CIRC-ENST0001.1-E2-E3')

    def test_to_convert_circ_rna_unknown_exon(self):
        """ circRNA with unkonwn exon. ExonNotFoundError should raise. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = self.create_base_record()
        record.start = 5
        record.end = 10
        record.exon_sizes = [5]
        record.exon_offsets = [0]
        with self.assertRaises(err.ExonNotFoundError):
            record.convert_to_circ_rna(anno)

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
        self.assertEqual(circ_record.id, 'CI-ENST0001.1-I1')

    def test_to_convert_ci_rna_negative(self):
        """ ciRNA """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)
        record = self.create_base_record()
        record.start = 12
        record.end = 17
        record.exon_sizes = [5]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        circ_record = record.convert_to_circ_rna(anno)
        self.assertIsInstance(circ_record, CircRNAModel)
        self.assertEqual(circ_record.id, 'CI-ENST0001.1-I2')

    def test_to_convert_ci_rna_fuzzy_end_positive(self):
        """ ciRNA with end before intron end positive strand """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = self.create_base_record()
        record.start = 12
        record.end = 15
        record.exon_sizes = [3]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        circ_record = record.convert_to_circ_rna(anno)
        self.assertIsInstance(circ_record, CircRNAModel)
        self.assertEqual(circ_record.id, 'CI-ENST0001.1-I1')

    def test_to_convert_ci_rna_fuzzy_end_negative(self):
        """ ciRNA with end before intron end negative strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)
        record = self.create_base_record()
        record.start = 14
        record.end = 17
        record.exon_sizes = [3]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        circ_record = record.convert_to_circ_rna(anno)
        self.assertIsInstance(circ_record, CircRNAModel)
        self.assertEqual(circ_record.id, 'CI-ENST0001.1-I2')

    def test_to_convert_ci_rna_invalid_end_positive(self):
        """ ciRNA with end after intron end on positive strand """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        record = self.create_base_record()
        record.start = 12
        record.end = 19
        record.exon_sizes = [7]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        with self.assertRaises(err.IntronNotFoundError):
            record.convert_to_circ_rna(anno)

    def test_to_convert_ci_rna_invalid_end_negative(self):
        """ ciRNA with end after intron end on negative strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(anno_data)
        record = self.create_base_record()
        record.start = 10
        record.end = 17
        record.exon_sizes = [7]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        with self.assertRaises(err.IntronNotFoundError):
            record.convert_to_circ_rna(anno)

    def test_to_convert_ci_rna_tx_unassociated_positive(self):
        """ The tx is not associated with the intron positive strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['transcripts'].append('ENST0002.1')
        tx2_attr = copy.deepcopy(ANNOTATION_ATTRS[1])
        tx2_attr['transcript_id'] = 'ENST0002.1'
        tx2_data = {
            'transcript_id': tx2_attr['transcript_id'],
            'chrom': 'chr1',
            'strand': 1,
            'transcript': (17, 35, tx2_attr),
            'exon': [(17, 23, tx2_attr), (27, 35, tx2_attr)]
        }
        anno_data['transcripts'].append(tx2_data)
        anno = create_genomic_annotation(anno_data)
        record = self.create_base_record()
        record.isoform_name = 'ENST0002.1'
        record.start = 12
        record.end = 17
        record.exon_sizes = [5]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        with self.assertRaises(err.IntronNotFoundError):
            record.convert_to_circ_rna(anno)

    def test_to_convert_ci_rna_tx_unassociated_negative(self):
        """ The tx is not associated with the intron on negative strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['transcripts'].append('ENST0002.1')
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        tx2_attr = copy.deepcopy(ANNOTATION_ATTRS[1])
        tx2_attr['transcript_id'] = 'ENST0002.1'
        tx2_data = {
            'transcript_id': tx2_attr['transcript_id'],
            'chrom': 'chr1',
            'strand': -1,
            'transcript': (5, 35, tx2_attr),
            'exon': [(5, 12, tx2_attr), (27, 35, tx2_attr)]
        }
        anno_data['transcripts'].append(tx2_data)
        anno = create_genomic_annotation(anno_data)
        record = self.create_base_record()
        record.isoform_name = 'ENST0002.1'
        record.start = 12
        record.end = 17
        record.exon_sizes = [5]
        record.exon_offsets = [0]
        record.circ_type = 'ciRNA'
        with self.assertRaises(err.IntronNotFoundError):
            record.convert_to_circ_rna(anno)
