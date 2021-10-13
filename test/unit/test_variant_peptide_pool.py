""" Test Module for VariantPeptidePool """
import unittest
from test.unit import create_aa_record, create_genomic_annotation
from test.unit.test_peptide_pool_splitter import ANNOTATION_DATA
from moPepGen.aa import VariantPeptidePool


class TestVariantPeptidePool(unittest.TestCase):
    """ Test cases for VariantPeptidePool """
    def test_add_peptide_label(self):
        """ Test the variant labels are kept when adding new peptide """
        data = [
            ['SSSSSSSSSR', 'ENST0001|SNV-100-A-T|1'],
            ['SSSSSSSSSR', 'ENST0002|SNV-200-G-C|17']
        ]
        peptides = [create_aa_record(*x) for x in data]
        pool = VariantPeptidePool(peptides={peptides[0]})
        pool.add_peptide(peptides[1], set())
        self.assertEqual(len(pool.peptides), 1)

        received = list(pool.peptides)[0].description
        expected = 'ENST0001|SNV-100-A-T|1 ENST0002|SNV-200-G-C|17'
        self.assertEqual(received, expected)

    def test_filter_base(self):
        """ test filter peptide pool """
        data = [
            ['SSSSSSSSSR', 'ENST0001|SNV-100-A-T|1'],
            ['SSSSSSSSAR', 'ENST0002|SNV-100-A-T|1'],
            ['SSSSSSSSCR', 'ENST0003|SNV-100-A-T|1'],
            ['SSSSSSSSGR', 'ENST0004|SNV-100-A-T|1'],
        ]
        exprs = {
            'ENST0001':10,
            'ENST0002':10,
            'ENST0003':5,
            'ENST0004':10,
        }
        anno = create_genomic_annotation(ANNOTATION_DATA)
        peptides = {create_aa_record(*x) for x in data}
        pool = VariantPeptidePool(peptides=peptides)
        filtered = pool.filter(
            exprs=exprs, cutoff=8, anno=anno, keep_all_coding=False,
            keep_all_noncoding=False
        )
        self.assertEqual(len(filtered.peptides), 3)
        seqs = {str(x.seq) for x in filtered.peptides}
        expected = {'SSSSSSSSSR', 'SSSSSSSSAR', 'SSSSSSSSGR'}
        self.assertEqual(seqs, expected)

    def test_filter_mixed(self):
        """ test filter when a peptide has two sources. """
        data = [
            ['SSSSSSSSSR', 'ENST0001|SNV-100-A-T|1'],
            ['SSSSSSSSAR', 'ENST0002|SNV-100-A-T|1'],
            ['SSSSSSSSCR', 'ENST0002|SNV-100-A-T|2 ENST0003|SNV-100-A-T|1'],
            ['SSSSSSSSGR', 'ENST0004|SNV-100-A-T|1'],
        ]
        exprs = {
            'ENST0001':10,
            'ENST0002':10,
            'ENST0003':5,
            'ENST0004':10,
        }
        anno = create_genomic_annotation(ANNOTATION_DATA)
        peptides = {create_aa_record(*x) for x in data}
        pool = VariantPeptidePool(peptides=peptides)
        filtered = pool.filter(
            exprs=exprs, cutoff=8, anno=anno, keep_all_coding=False,
            keep_all_noncoding=False
        )
        self.assertEqual(len(filtered.peptides), 4)

        exprs['ENST0002'] = 4
        filtered = pool.filter(
            exprs=exprs, cutoff=8, anno=anno, keep_all_coding=False,
            keep_all_noncoding=False
        )
        self.assertEqual(len(filtered.peptides), 2)

    def test_filter_keep_all_coding(self):
        """ test filter peptide pool keep all coding """
        data = [
            ['SSSSSSSSSR', 'ENST0001|SNV-100-A-T|1'],
            ['SSSSSSSSAR', 'ENST0002|SNV-100-A-T|1'],
            ['SSSSSSSSCR', 'ENST0003|SNV-100-A-T|1'],
            ['SSSSSSSSGR', 'ENST0004|SNV-100-A-T|1'],
        ]
        exprs = {
            'ENST0001':10,
            'ENST0002':10,
            'ENST0003':5,
            'ENST0004':10,
        }
        anno = create_genomic_annotation(ANNOTATION_DATA)
        peptides = {create_aa_record(*x) for x in data}
        pool = VariantPeptidePool(peptides=peptides)
        filtered = pool.filter(
            exprs=exprs, cutoff=8, anno=anno, keep_all_coding=True,
            keep_all_noncoding=False
        )
        self.assertEqual(len(filtered.peptides), 4)

    def test_filter_keep_all_noncoding(self):
        """ test filter peptide pool keep all coding """
        data = [
            ['SSSSSSSSSR', 'ENST0001|SNV-100-A-T|1'],
            ['SSSSSSSSAR', 'ENST0002|SNV-100-A-T|1'],
            ['SSSSSSSSCR', 'ENST0003|SNV-100-A-T|1'],
            ['SSSSSSSSGR', 'ENST0004|SNV-100-A-T|1'],
        ]
        exprs = {
            'ENST0001':10,
            'ENST0002':10,
            'ENST0003':10,
            'ENST0004':5,
        }
        anno = create_genomic_annotation(ANNOTATION_DATA)
        peptides = {create_aa_record(*x) for x in data}
        pool = VariantPeptidePool(peptides=peptides)
        filtered = pool.filter(
            exprs=exprs, cutoff=8, anno=anno, keep_all_coding=False,
            keep_all_noncoding=True
        )
        self.assertEqual(len(filtered.peptides), 4)
