""" Test the CLI for callNoncoding """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestCallNoncodingPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callNoncoding """

    def test_call_noncoding_peptides(self):
        """ test call noncoding peptides """
        args = argparse.Namespace()
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.inclusion_biotypes = None
        args.exclusion_biotypes = None
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_noncoding_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)