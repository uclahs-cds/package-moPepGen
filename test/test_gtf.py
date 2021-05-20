""" Test the GTF files are loaded and handled properly
"""
from unittest import TestCase
from moPepGen.gtf import TranscriptGTFDict
from moPepGen.dna.DNASeqDict import DNASeqDict


class TestGTF(TestCase):
    """ Test case for GTF modules
    """
    def load_gtf(self, path:str):
        gtf = TranscriptGTFDict()
        gtf.dump_gtf(path)
        return gtf
    
    def test_dump_gtf(self):
        """ Test that the VapIO.parse returns an iterabale of vep_record
        """
        gtf = self.load_gtf('test/files/gtf_example_gencode.gtf')
        
        self.assertIsInstance(gtf, TranscriptGTFDict)
        self.assertEqual(len(gtf), 5)
        for key, val in gtf.items():
            self.assertEqual(val.transcript.attributes['transcript_id'], key)
            for cds in val.cds:
                self.assertEqual(cds.attributes['transcript_id'], key)
            for exon in val.exon:
                self.assertEqual(exon.attributes['transcript_id'], key)
    

class TestTranscriptAnnotationDict(TestCase):
    """ Test case for Transcript Annotation Data Models """
    def test_get_cds_sequence(self):
        genome = DNASeqDict()
        genome.dump_fasta(
            'test/files/downsampled_set/gencode_v34_genome_chr22.fasta')
        gtf = TranscriptGTFDict()
        gtf.dump_gtf('test/files/downsampled_set/gencode_v34_chr22.gtf')

        cdna_seq = gtf['ENST00000543038.1'].get_cdna_sequence(genome['chr22'])
        self.assertEqual(len(cdna_seq.seq), len(cdna_seq.locations[0].query))
        self.assertEqual(len(cdna_seq.seq), len(cdna_seq.locations[0].ref))