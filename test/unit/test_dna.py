""" Test the modules that handles the DNA sequences """
import unittest
from test.unit import create_dna_seq_with_coordinates
from Bio.Seq import Seq
from moPepGen import dna
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation


SEQ = 'CCCTACTGGTCCTTCTGCCTTAGCCACAGGTTCTGAAACCAAAGCAAAACCACCAGAGAG' +\
'TGATTCATTCCCGTGGAGTCTCCATCTGAGCCCTTTCCTAGTCCAGGCATCCCGATGTTG' +\
'GTGGATGGCCCATCTGAGCGGCCAGCCCTGTGCTTCTTGCTGTTGGCTGTGGCAATGTCT' +\
'TTCTTC'

def load_fasta_from_disk(path:str):
    """ Load a fasta file from disk """
    seqs = dna.DNASeqDict()
    seqs.dump_fasta(path)
    return seqs

class TestDNASeqRecord(unittest.TestCase):
    """ Test DNA """
    def test_dump_fasta(self):
        """ Test that the fasta file can be loaded """
        seqs = load_fasta_from_disk('test/files/genome.fasta')
        self.assertIsInstance(seqs, dna.DNASeqDict)
        for record in seqs.values():
            self.assertIsInstance(record, dna.DNASeqRecord)

    def test_dna_seq_dict_type(self):
        """ Test that only dna.DNASeqRecord is allowed. """
        seqs = load_fasta_from_disk('test/files/genome.fasta')

        with self.assertRaises(TypeError):
            seqs[list(seqs.keys())[0]] = 1

    def test_find_orf(self):
        """ Test searching ORF """
        seq = dna.DNASeqRecord(Seq(SEQ))
        self.assertEqual(seq.find_orf_first(), 114)
        for i in seq.find_orf_all():
            self.assertEqual(seq.seq[i:i+3], 'ATG')

    def test_add_dna_record(self):
        """ test DNASeq.__add__ """
        seq1 = dna.DNASeqRecord(Seq('CCCTACTGGT'))
        seq2 = dna.DNASeqRecord(Seq('CCTTCTGCCT'))
        seq = seq1 + seq2
        self.assertIsInstance(seq, dna.DNASeqRecord)
        self.assertEqual(str(seq.seq), 'CCCTACTGGTCCTTCTGCCT')

    def test_find_all_start_codons(self):
        """ find all start codons """
        seq = dna.DNASeqRecord('CCCCCCCCATGCCCCCCCCATGCCCCCCCATGCCCCCCCC')
        starts = seq.find_all_start_codons()
        self.assertEqual(len(starts), 3)
        self.assertTrue(all(seq.seq[i:i+3] == 'ATG' for i in starts))

class TestDNASeqWithCoordinates(unittest.TestCase):
    """ Test case for moPepGen.dna.DNASeqRecordWithCoordinates """
    def test_get_item_case(self):
        """ test __getitem__ """

        location = MatchedLocation(
            query=FeatureLocation(start=0, end=20),
            ref=FeatureLocation(start=0, end=20)
        )
        seq = dna.DNASeqRecordWithCoordinates(
            seq=Seq('CCCTACTGGTCCTTCTGCCT'),
            locations=[location]
        )
        seq1 = seq[:10]
        self.assertEqual(str(seq1.seq), 'CCCTACTGGT')
        self.assertEqual(int(seq1.locations[0].ref.start), 0)
        self.assertEqual(int(seq1.locations[0].ref.end), 10)
        self.assertEqual(int(seq1.locations[0].query.start), 0)
        self.assertEqual(int(seq1.locations[0].query.end), 10)

        seq1 = seq[10:19]
        self.assertEqual(str(seq1.seq), 'CCTTCTGCC')
        self.assertEqual(int(seq1.locations[0].ref.start), 10)
        self.assertEqual(int(seq1.locations[0].ref.end), 19)
        self.assertEqual(int(seq1.locations[0].query.start), 0)
        self.assertEqual(int(seq1.locations[0].query.end), 9)

    def test_add(self):
        """ test __add__ """
        loc1 = [((0,10), (0,10))]
        loc2 = [((0,10), (20,30))]

        seq1 = create_dna_seq_with_coordinates('CCCTACTGGT', loc1, [0, None])
        seq2 = create_dna_seq_with_coordinates('CCTTCTGCCT', loc2, [0, None])

        seq = seq1 + seq2
        self.assertEqual(str(seq.seq), 'CCCTACTGGTCCTTCTGCCT')
        self.assertEqual(int(seq.locations[0].ref.start), 0)
        self.assertEqual(int(seq.locations[0].ref.end), 10)
        self.assertEqual(int(seq.locations[1].ref.start), 20)
        self.assertEqual(int(seq.locations[1].ref.end), 30)
        self.assertEqual(int(seq.locations[0].query.start), 0)
        self.assertEqual(int(seq.locations[0].query.end), 10)
        self.assertEqual(int(seq.locations[1].query.start), 10)
        self.assertEqual(int(seq.locations[1].query.end), 20)

        loc1 = [((0,10), (0,10))]
        loc2 = [((1,10), (10,19))]

        seq1 = create_dna_seq_with_coordinates('CCCTACTGGT', loc1, [0, None])
        seq2 = create_dna_seq_with_coordinates('CCTTCTGCCT', loc2, [0, None])
        seq = seq1 + seq2
        self.assertEqual(len(seq.locations), 2)
        self.assertEqual(int(seq.locations[1].query.start), 11)
        self.assertEqual(int(seq.locations[1].query.end), 20)

        loc1 = [((0,10), (0,10))]
        loc2 = [((0,9), (11,20))]

        seq1 = create_dna_seq_with_coordinates('CCCTACTGGT', loc1, [0, None])
        seq2 = create_dna_seq_with_coordinates('CCTTCTGCCT', loc2, [0, None])
        seq = seq1 + seq2
        self.assertEqual(len(seq.locations), 2)
        self.assertEqual(int(seq.locations[1].query.start), 10)
        self.assertEqual(int(seq.locations[1].query.end), 19)

    def test_getitem(self):
        """ test get item """
        loc1 = [((0,4),(0,4))]
        seq1 = create_dna_seq_with_coordinates('AAAC', loc1, (0,None))
        seq = seq1[1:3]
        self.assertEqual(str(seq.seq), 'AA')
        self.assertEqual(len(seq.locations), 1)
        self.assertEqual(seq.locations[0].query.start, 0)
        self.assertEqual(seq.locations[0].query.end, 2)
        self.assertEqual(seq.locations[0].ref.start, 1)
        self.assertEqual(seq.locations[0].ref.end, 3)

        loc1 = [((0,4),(0,4)), ((5,8), (5,8))]
        seq1 = create_dna_seq_with_coordinates('AAACGGGT', loc1, (0,None))
        seq = seq1[3:7]
        self.assertEqual(len(seq.locations), 2)
        self.assertEqual(seq.locations[0].query.start, 0)
        self.assertEqual(seq.locations[0].query.end, 1)
        self.assertEqual(seq.locations[0].ref.start, 3)
        self.assertEqual(seq.locations[0].ref.end, 4)
        self.assertEqual(seq.locations[1].query.start, 2)
        self.assertEqual(seq.locations[1].query.end, 4)
        self.assertEqual(seq.locations[1].ref.start, 5)
        self.assertEqual(seq.locations[1].ref.end, 7)

        seq = seq1[4:7]
        self.assertEqual(len(seq.locations), 1)
        self.assertEqual(seq.locations[0].query.start, 1)
        self.assertEqual(seq.locations[0].query.end, 3)
        self.assertEqual(seq.locations[0].ref.start, 5)
        self.assertEqual(seq.locations[0].ref.end, 7)

        seq = seq1[2:5]
        self.assertEqual(len(seq.locations), 1)
        self.assertEqual(seq.locations[0].query.start, 0)
        self.assertEqual(seq.locations[0].query.end, 2)
        self.assertEqual(seq.locations[0].ref.start, 2)
        self.assertEqual(seq.locations[0].ref.end, 4)

    def test_get_query_index(self):
        """ Test the correct query index is returned. """
        location = MatchedLocation(
            query=FeatureLocation(start=0, end=20),
            ref=FeatureLocation(start=101, end=121)
        )
        seq = dna.DNASeqRecordWithCoordinates(
            seq='ACCGTTGAGTTGGGTACCGC',
            locations=[location]
        )
        self.assertEqual(seq.get_query_index(105), 4)
