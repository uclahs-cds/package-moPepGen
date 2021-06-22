""" Module for testing seqvar.io """
import unittest
from moPepGen import seqvar


class TestSeqvarIO(unittest.TestCase):
    """ Test case for parsing seqvar files. """
    def test_seqvar_parse(self):
        """ Test parsing seqvar files. """
        mop_path = 'test/files/vep_test_files/vep_moPepGen.txt'
        i = 0
        for record in seqvar.io.parse(mop_path):
            i += 1
            self.assertIsInstance(record, seqvar.VariantRecord)
            if i > 5:
                break
