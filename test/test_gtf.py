""" Test the GTF files are loaded and handled properly
"""
from unittest import TestCase
from moPepGen.TranscriptGTFDict import TranscriptGTFDict
from moPepGen.GTFRecord import GTFRecord


class TestGTF(TestCase):
    """ Test case for GTF modules
    """
    def load_gtf_from_disk(self, path:str):
        gtf = TranscriptGTFDict()
        gtf.dump_gtf(path)
        return gtf
    
    def test_dump_gtf(self):
        """ Test that the VapIO.parse returns an iterabale of vep_record
        """
        gtf = self.load_gtf_from_disk('test/files/gtf_example_gencode.gtf')
        
        self.assertIsInstance(gtf, TranscriptGTFDict)
        self.assertEqual(len(gtf), 6)
        for key, val in gtf.items():
            self.assertEqual(val.transcript.attributes['transcript_id'], key)
            for cds in val.cds:
                self.assertEqual(cds.attributes['transcript_id'], key)
            for exon in val.exon:
                self.assertEqual(exon.attributes['transcript_id'], key)
    
