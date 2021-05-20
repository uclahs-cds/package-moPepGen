""" Test the modules that handles the DNA sequences """
from moPepGen.dna.DNASeqRecord import DNASeqRecordWithCoordinates
from unittest import TestCase
from moPepGen.DNASeqDict import DNASeqDict, DNASeqRecord


SEQ = 'CCCTACTGGTCCTTCTGCCTTAGCCACAGGTTCTGAAACCAAAGCAAAACCACCAGAGAG' +\
'TGATTCATTCCCGTGGAGTCTCCATCTGAGCCCTTTCCTAGTCCAGGCATCCCGATGTTG' +\
'GTGGATGGCCCATCTGAGCGGCCAGCCCTGTGCTTCTTGCTGTTGGCTGTGGCAATGTCT' +\
'TTCTTC'

class TestDNASeq(TestCase):
    """ Test DNA """
    def load_fasta_from_disk(self, path:str):
        dna = DNASeqDict()
        dna.dump_fasta(path)
        return dna
    
    def test_dump_fasta(self):
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
        
    def test_find_orf(self):
        """ Test searching ORF """
        seq = DNASeqRecord(seq)
        self.assertEqual(seq.find_orf_first(), 114)
        for i in seq.find_orf_all():
            self.assertEqual(seq.seq[i:i+3], 'ATG')
            