""" Test the command line interface """
import argparse
from pathlib import Path
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseVEP(TestCaseIntegration):
    """ Test cases for moPepGen parseVEP """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base args """
        args = argparse.Namespace()
        args.command = 'parseVEP'
        args.input_path = []
        args.index_dir = None
        args.source = 'gSNP'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.output_path = Path(self.work_dir/'vep.gvf')
        args.quiet = True
        return args

    def test_parse_vep(self):
        """ Test parsing VEP output into GVF """
        args = self.create_base_args()
        args.input_path = [
            self.data_dir/'vep'/'vep_snp.txt',
            self.data_dir/'vep'/'vep_indel.txt'
        ]
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_vep_gz(self):
        """ Test parsing gzipped VEP output into GVF """
        args = self.create_base_args()
        args.input_path = [
            self.data_dir/'vep'/'vep_snp.txt.gz',
            self.data_dir/'vep'/'vep_indel.txt.gz'
        ]
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)
