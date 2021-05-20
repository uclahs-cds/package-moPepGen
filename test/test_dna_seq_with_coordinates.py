""" Test DNASeqWithLocation """
from unittest import TestCase
from moPepGen.dna import DNASeqRecordWithCoordinates, DNASeqRecord, \
    MatchedLocation
from moPepGen.dna.DNASeqDict import DNASeqDict
from moPepGen.gtf.TranscriptGTFDict import TranscriptGTFDict
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

    def test_get_item(self):
        """ Test """
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
    
    def test_enzymatic_cleave(self):
        genome = DNASeqDict()
        genome.dump_fasta(
            'test/files/downsampled_set/gencode_v34_genome_chr22.fasta')
        gtf = TranscriptGTFDict()
        gtf.dump_gtf('test/files/downsampled_set/gencode_v34_chr22.gtf')

        cdna_seq = gtf['ENST00000543038.1'].get_cdna_sequence(genome['chr22'])
        cdna_seq.enzymatic_cleave(rule='trypsin')

    def test_get_query_index(self):
        location = MatchedLocation(
            query=FeatureLocation(start=0, end=20),
            ref=FeatureLocation(start=101, end=121)
        )
        seq = DNASeqRecordWithCoordinates(
            seq='ACCGTTGAGTTGGGTACCGC',
            locations=[location]
        )
        self.assertEqual(seq.get_query_index(105), 4)