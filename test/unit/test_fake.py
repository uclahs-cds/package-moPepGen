""" Unit test for the fake module """
import unittest
import random
from moPepGen import fake


class TestCaseFake(unittest.TestCase):
    """ Test cases for fake """
    def test_fake_transcript_model(self):
        """ Test fake_transcript_model """
        tx_model = fake.fake_transcript_model(
            n_exons=2, is_coding=True, chrom='chrF', strand=1, start_pos=0,
            gene_id='FAKEG000001', transcript_id='FAKET000001',
            protein_id='FAKEP000001', cds_start_nf=False, mrna_end_nf=False,
            min_exon_size=10, max_exon_size=300, min_intron_size=20,
            max_intron_size=500
        )
        self.assertEqual(len(tx_model.exon), 2)
        self.assertEqual(tx_model.transcript.strand, 1)
        self.assertTrue(all(10 <= len(x) <= 300 for x in tx_model.exon))

    def test_fake_genome_and_annotation(self):
        """ Test fake_genome_and_annotation """
        random.seed(123124)
        anno = fake.fake_genomic_annotation(
            n_genes=4, chrom='chrF', min_exons=1, max_exons=10, min_exon_size=10,
            max_exon_size=300, min_intron_size=20, max_intron_size=500,
            min_intergenic_size=10, max_intergenic_size=50
        )
        self.assertEqual(len(anno.genes), 4)

        genome = fake.fake_genome(anno)
        self.assertEqual(len(genome), 1)

        for tx_model in anno.transcripts.values():
            if tx_model.is_protein_coding:
                chrom = tx_model.transcript.chrom
                tx_seq = tx_model.get_transcript_sequence(genome[chrom])
                aa_seq = tx_seq.seq[tx_seq.orf.start:tx_seq.orf.end].translate()
                if not tx_model.is_cds_start_nf():
                    self.assertTrue(aa_seq.startswith('M'))
                self.assertTrue('*' not in aa_seq)
