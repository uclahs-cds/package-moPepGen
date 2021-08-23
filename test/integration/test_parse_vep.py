""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseVEP(TestCaseIntegration):
    """ Test cases for moPepGen parseVEP """

    def test_parse_vep(self):
        """ Test genreate index """
        args = argparse.Namespace()
        args.vep_txt = [
            self.data_dir/'vep'/'vep_snp.txt',
            self.data_dir/'vep'/'vep_indel.txt'
        ]
        # args.index_dir = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'vep')
        args.verbose = False
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep.gvf'}
        self.assertEqual(files, expected)
