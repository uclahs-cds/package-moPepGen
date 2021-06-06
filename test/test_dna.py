""" Test the modules that handles the DNA sequences """
from unittest import TestCase
from Bio.Seq import Seq
from moPepGen.dna import DNASeqDict, DNASeqRecord


SEQ = 'CCCTACTGGTCCTTCTGCCTTAGCCACAGGTTCTGAAACCAAAGCAAAACCACCAGAGAG' +\
'TGATTCATTCCCGTGGAGTCTCCATCTGAGCCCTTTCCTAGTCCAGGCATCCCGATGTTG' +\
'GTGGATGGCCCATCTGAGCGGCCAGCCCTGTGCTTCTTGCTGTTGGCTGTGGCAATGTCT' +\
'TTCTTC'

def load_fasta_from_disk(path:str):
    """ Load a fasta file from disk """
    dna = DNASeqDict()
    dna.dump_fasta(path)
    return dna

class TestDNASeq(TestCase):
    """ Test DNA """
    def test_dump_fasta(self):
        """ Test that the fasta file can be loaded """
        dna = load_fasta_from_disk('test/files/genome_example.fa')
        self.assertIsInstance(dna, DNASeqDict)
        for record in dna.values():
            self.assertIsInstance(record, DNASeqRecord)

    def test_dna_seq_dict_type(self):
        """ Test that only DNASeqRecord is allowed. """
        dna = load_fasta_from_disk('test/files/genome_example.fa')

        with self.assertRaises(TypeError):
            dna[list(dna.keys())[0]] = 1

    def test_find_orf(self):
        """ Test searching ORF """
        seq = DNASeqRecord(Seq(SEQ))
        self.assertEqual(seq.find_orf_first(), 114)
        for i in seq.find_orf_all():
            self.assertEqual(seq.seq[i:i+3], 'ATG')
