""" Module to test CircRNAModel """
import unittest
from moPepGen import CircRNA


class TestCircRNA(unittest.TestCase):
    """ Test case for circ RNA """
    def test_parse_circ_rna_bed(self):
        """ Test to parse circRNA bed file """
        circs = list(CircRNA.parse('test/files/circRNA/circ_rna.tsv'))
        for circ in circs:
            self.assertIsInstance(circ, CircRNA.CircRNAModel)
        self.assertEqual(len(circs[0].fragments), 2)
