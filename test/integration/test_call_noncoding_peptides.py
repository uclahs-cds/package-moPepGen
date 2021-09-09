""" Test the CLI for callNoncoding """
import argparse
import os
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli


def create_base_args() -> argparse.Namespace:
    """ Create a base args """
    args = argparse.Namespace()
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.output_fasta = None
    args.inclusion_biotypes = None
    args.exclusion_biotypes = None
    args.min_tx_length = 21
    args.cleavage_rule = 'trypsin'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.verbose = False
    return args

class TestCallNoncodingPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callNoncoding """

    def test_call_noncoding_peptides_case1(self):
        """ test call noncoding peptides """
        args = create_base_args()
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_fasta = self.work_dir/'noncoding.fasta'
        cli.call_noncoding_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding.fasta'}
        self.assertEqual(files, expected)
        with open(self.work_dir/'noncoding.fasta', 'rt') as handle:
            peptides = list(SeqIO.parse(handle, 'fasta'))
            ids = [p.id for p in peptides]
            self.assertEqual(len(ids),len(set(ids)))

    def test_call_noncoding_peptides_case2(self):
        """ test call noncoding peptides when no ORF is found """
        args = create_base_args()
        ref_dir = self.data_dir/'downsampled_reference/ENST00000644482.1'
        args.genome_fasta = ref_dir/'genome.fasta'
        args.annotation_gtf = ref_dir/'annotation.gtf'
        args.proteome_fasta = ref_dir/'proteome.fasta'
        args.output_fasta = self.work_dir/'noncoding.fasta'
        cli.call_noncoding_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding.fasta'}
        self.assertEqual(files, expected)
        size = os.stat(self.work_dir/'noncoding.fasta').st_size
        self.assertEqual(size, 0)
