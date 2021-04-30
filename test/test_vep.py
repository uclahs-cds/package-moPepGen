""" Test the moPepGen.io.VepIO module
"""
from unittest import TestCase
from moPepGen.TranscriptVEPDict import TranscriptVEPDict


class TestVep(TestCase):
    """ Test case for VEP modules
    """

    def load_vep_from_disk(self, file: str):
        """ Load the test VEP file from disk """
        vep_path = file
        vep = TranscriptVEPDict()
        vep.dump_vep(vep_path)
        return vep

    def test_vep_from_file(self):
        """ Test that the VapIO.parse returns an iterabale of vep_record
        """
        vep = self.load_vep_from_disk('test/files/vep_example_snp.txt')
        self.assertIsInstance(vep, TranscriptVEPDict)
        for key, val in vep.items():
            for vep_record in val:
                self.assertEqual(vep_record.feature, key)
    
    def test_transcript_vep_records_type(self):
        """ Check that TYpeError is raised then trying to assign elements with
        wrong type.
        """
        vep = self.load_vep_from_disk('test/files/vep_example_snp.txt')
        vep_records = list(vep.values())[0]

        with self.assertRaises(TypeError):
            vep_records[0] = 1
        
        with self.assertRaises(TypeError):
            vep_records.append('string')
        
        with self.assertRaises(TypeError):
            vep_records.insert(1, 1.25)
        
        with self.assertRaises(TypeError):
            vep_records.extend(1, [1,2,3])
