""" Test DNASeqWithLocation """
from unittest import TestCase
from moPepGen.dna import DNASeqRecordWithCoordinates, \
    MatchedLocation
from moPepGen.SeqFeature import FeatureLocation


SEQ = 'ATGTTGGTGGATGGCCCATCTGAGCGGCCAGCCCTGTGCTTCTTGCTGTTGGCTGTGGCAAT' +\
    'GTCTTTCTTCGGCTCAGCTCTATCCATAGATGAAACACGGGCGCATCTGTTGTTGAAAGAAAAG' +\
    'ATGATGCGGCTGGGGGGGCGGCTGGTGCTGAACACCAAGGAGGAGCTGGCCAATGAGAGGCTCA' +\
    'TGACGCTCAAAATCGCTGAGATGAAGGAGGCCATGAGGACCCTGATATTCCCACCCAGCATGCA' +\
    'CTTTTTCCAGGCCAAGCATCTCATTGAGAGAAGTCAAGTGTTTAATATTCTAAGGATGATGCCA' +\
    'AAAGGGGCTGCCTTGCACCTCCATGACATTGGCATCGTGACTATGGACTGGCTGGTGAGGAATG' +\
    'TCACCTACAGGCCTCACTGCCACATCTGTTTCACCCCAAGGGGGATCATGCAGTTCAGATTTGCT' +\
    'CACCCAACT'
TRANCRIPT_ID = 'ENST00000543038.1'


class TestDNASeqRecordWithCoordinates(TestCase):
    """ Test case for DNARecordWithCoordinate """
    def test_get_item(self):
        """ Test __getitem__ """
        location = MatchedLocation(
            query=FeatureLocation(start=0, end=20),
            ref=FeatureLocation(start=101, end=121)
        )
        seq = DNASeqRecordWithCoordinates(
            seq='ACCGTTGAGTTGGGTACCGC',
            locations=[location]
        )
        subseq = seq[5:15]
        self.assertEqual(str(subseq.seq), 'TGAGTTGGGT')
        self.assertEqual(len(subseq.locations), 1)
        self.assertEqual(subseq.locations[0].ref.start, 106)
        self.assertEqual(subseq.locations[0].ref.end, 116)

    def test_get_query_index(self):
        """ Test the correct query index is returned. """
        location = MatchedLocation(
            query=FeatureLocation(start=0, end=20),
            ref=FeatureLocation(start=101, end=121)
        )
        seq = DNASeqRecordWithCoordinates(
            seq='ACCGTTGAGTTGGGTACCGC',
            locations=[location]
        )
        self.assertEqual(seq.get_query_index(105), 4)
