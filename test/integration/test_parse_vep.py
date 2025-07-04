""" Test the command line interface """
import argparse
from pathlib import Path
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseVEP(TestCaseIntegration):
    """ Test cases for moPepGen parseVEP """
    def test_parse_vep_cli(self):
        """ test parseVEP cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli parseVEP \\
            -i {self.data_dir}/vep/vep_snp.txt \\
            -o {self.work_dir}/vep.gvf \\
            -g {self.data_dir}/genome.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta \\
            --source VEP
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

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
        args.reference_source = None
        args.output_path = Path(self.work_dir/'vep.gvf')
        args.output_prefix = Path(self.work_dir/'vep')
        args.samples = []
        args.quiet = True
        return args

    def test_parse_vep(self):
        """ Test parsing VEP output into GVF """
        args = self.create_base_args()
        args.input_path = [
            #self.data_dir/'vep'/'vep_snp.txt',
            self.data_dir/'vep'/'vep_indel.txt'
        ]
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_vep2(self):
        """ Failed records are skipped """
        args = self.create_base_args()
        args.input_path = [
            self.data_dir/'vep/vep_indel2.txt'
        ]
        args.skip_failed = True
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

    def test_parse_vep_vcf(self):
        """ Test parsing VEP output in VCF format into GVF """
        args = self.create_base_args()
        args.input_path = [
            self.data_dir/'vep'/'vep_snp.vcf'
        ]
        cli.parse_vep(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_UCLA0001.gvf'}
        self.assertEqual(files, expected)
        for file_name in expected:
            self.assert_gvf_order(self.work_dir/file_name, args.annotation_gtf)
