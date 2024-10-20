""" Integration test for parseCIRCexplorer """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseCIRCexplorer(TestCaseIntegration):
    """ Test cases for moPepGen parseCIRCexplorer """
    def test_parse_circexplorer2(self):
        """ Test parseCIRCexplorer """
        args = argparse.Namespace()
        args.command = 'parseCIRCexplorer'
        args.input_path = self.data_dir/'circRNA/CIRCexplorer_circularRNA_known.txt'
        args.output_path = self.work_dir/'circ.gvf'
        args.source = 'circRNA'
        args.index_dir = None
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.reference_source = None
        args.circexplorer3 = False
        args.min_read_number = 1
        args.intron_start_range = '-2,0'
        args.intron_end_range = '-100,2'
        args.quiet = True
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(self.work_dir/'circ.gvf', args.annotation_gtf)

    def test_parse_circexplorer3(self):
        """ Test parseCIRCexplorer for CIRCexplorer3 """
        args = argparse.Namespace()
        args.command = 'parseCIRCexplorer'
        args.input_path = self.data_dir/'circRNA/CIRCexplorer3_circularRNA_known.txt'
        args.output_path = self.work_dir/'circ.gvf'
        args.source = 'circRNA'
        args.index_dir = None
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.reference_source = None
        args.circexplorer3 = True
        args.min_read_number = 1
        args.min_fbr_circ = None
        args.min_circ_score = None
        args.intron_start_range = '-2,0'
        args.intron_end_range = '-100,2'
        args.quiet = True
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

        args.min_fbr_circ = 1
        args.min_circ_score = 1
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_circexplorer_cli(self):
        """ Test parseCIRCexplorer cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli parseCIRCexplorer \\
            -i {self.data_dir}/circRNA/CIRCexplorer3_circularRNA_known.txt \\
            -o {self.work_dir}/circ.gvf \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta \\
            --source circRNA
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise
