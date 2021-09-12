""" Test Module for VariantPeptidePool """
import unittest
from test.unit import create_aa_record
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
