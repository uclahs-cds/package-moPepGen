""" Module to test CircRNAModel """
import unittest
from moPepGen import circ


class TestCircRNA(unittest.TestCase):
    """ Test case for circ RNA """
    def test_parse_circ_rna_bed(self):
        """ Test to parse circRNA bed file """
        records = list(circ.io.parse('test/files/circRNA/circ_rna.tsv'))
        for record in records:
            self.assertIsInstance(record, circ.io.CircRNAModel)
        self.assertEqual(len(records[0].fragments), 2)
