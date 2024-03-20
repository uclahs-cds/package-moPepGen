""" Test the CLI for callAltStart """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


def create_base_args() -> argparse.Namespace:
    """ Create a base args """
    args = argparse.Namespace()
    args.command = 'callAltTranslation'
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.reference_source = None
    args.output_path = None
    args.output_orf = None
    args.w2f_reassignment = True
    args.orf_assignment = 'max'
    args.cleavage_rule = 'trypsin'
    args.cleavage_exception = 'trypsin_exception'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.quiet = True
    return args

class TestCallAltStart(TestCaseIntegration):
    """ Test cases for moPepGen callAltTranslation """

    def test_call_alt_start_cli(self):
        """ test callAltTranslation cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli callAltStart \\
            -o {self.work_dir}/alt_start_peptide.fasta \\
            --output-orf {self.work_dir}/alt_start_orf.fasta \\
            -g {self.data_dir}/genome.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_call_alt_start_case1(self):
        """ test call alt translation peptides """
        args = create_base_args()
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'alt_start_peptides.fasta'
        args.output_orf = self.work_dir/'alt_start_orf.fasta'

        cli.call_alt_start(args)

        self.assertTrue(args.output_path.exists())
        self.assertTrue(args.output_orf.exists())
