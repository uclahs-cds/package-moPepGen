""" Test the command line interface """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestGenerateIndex(TestCaseIntegration):
    """ Test cases for moPepGen generateIndex """
    def test_generate_index_cli(self):
        """ Test generateIndex cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli generateIndex \\
            -o {self.work_dir}/index \\
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

    def test_generate_index_case1(self):
        """ Test genreate index """
        args = argparse.Namespace()
        args.genome_fasta = self.data_dir / 'genome.fasta'
        args.annotation_gtf = self.data_dir / 'annotation.gtf'
        args.proteome_fasta = self.data_dir / 'translate.fasta'
        args.invalid_protein_as_noncoding = False
        args.cleavage_rule = 'trypsin'
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.miscleavage = 2
        args.quiet = True
        args.output_dir = self.work_dir / 'index'
        args.output_dir.mkdir(parents=False, exist_ok=True)
        cli.generate_index(args)
        files = {str(file.name) for file in args.output_dir.glob('*.pkl')}
        expected = {'genome.pkl', 'proteome.pkl', 'annotation.pkl',
            'canonical_peptides.pkl', 'coding_transcripts.pkl'}
        self.assertEqual(files, expected)
