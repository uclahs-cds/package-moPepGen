""" Test module for parseREDItools """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseREDItools(TestCaseIntegration):
    """ Test parseREDItools """
    def test_parse_reditools_cli(self):
        """ Test parseREDItools cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli parseREDItools \\
            -i {self.data_dir}/reditools/reditools_annotated.txt \\
            -o {self.work_dir}/reditools.gvf \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta \\
            --source RES
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_parse_reditools_case1(self):
        """ Test parse reditools """
        args = argparse.Namespace()
        args.command = 'parseREDItools'
        args.source = 'RNAEditingSite'
        args.input_path = self.data_dir/'reditools/reditools_annotated.txt'
        args.transcript_id_column = 17
        args.index_dir = None
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.reference_source = None
        args.output_path = self.work_dir/'reditools.gvf'
        args.quiet = True
        args.min_coverage_alt = 3
        args.min_frequency_alt = 0.1
        args.min_coverage_rna = 10
        args.min_coverage_dna = 10
        cli.parse_reditools(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'reditools.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(self.work_dir/'reditools.gvf', args.annotation_gtf)
