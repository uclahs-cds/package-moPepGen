""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseVEP(TestCaseIntegration):
    """ Test cases for moPepGen parseVEP """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base args """
        args = argparse.Namespace()
        args.command = 'parseVEP'
        args.vep_txt = []
        args.index_dir = None
        args.source = 'gSNP'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.output_prefix = str(self.work_dir/'vep')
        args.quiet = True
        return args

    def test_parse_vep(self):
        """ Test parsing VEP output into GVF """
        args = self.create_base_args()
        args.vep_txt = [
            self.data_dir/'vep'/'vep_snp.txt',
            self.data_dir/'vep'/'vep_indel.txt'
        ]
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep.gvf'}
        self.assertEqual(files, expected)

    def test_parse_vep_gz(self):
        """ Test parsing gzipped VEP output into GVF """
        args = self.create_base_args()
        args.vep_txt = [
            self.data_dir/'vep'/'vep_snp.txt.gz',
            self.data_dir/'vep'/'vep_indel.txt.gz'
        ]
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep.gvf'}
        self.assertEqual(files, expected)
