""" Test the modules that handles the DNA sequences """
from unittest import TestCase
from moPepGen.DNASeqDict import DNASeqDict, DNASeqRecord


class TestDNASeq(TestCase):
    """ Test DNA """
    def load_fasta_from_disk(self, path:str):
        dna = DNASeqDict()
        dna.dump_genome(path)
        return dna
    
    def test_dump_genome(self):
        """ Test that the fasta file can be loaded """
        dna = self.load_fasta_from_disk('test/files/genome_example.fa')
        self.assertIsInstance(dna, DNASeqDict)
        for record in dna.values():
            self.assertIsInstance(record, DNASeqRecord)
    
    def test_dna_seq_dict_type(self):
        """ Test that only DNASeqRecord is allowed. """
        dna = self.load_fasta_from_disk('test/files/genome_example.fa')

        with self.assertRaises(TypeError):
            dna[list(dna.keys())[0]] = 1
        