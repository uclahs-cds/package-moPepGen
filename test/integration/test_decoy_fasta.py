""" Integration test for decoyFasta """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration, TestFastaWriterMixin
from Bio import SeqIO
from moPepGen import cli


TARGET_DB = [
    ('MIPTGGERQWQLIPPVK', 'ENST0001|SNV-100-A-T'),
    ('IRPPGGQIICCVCSSSR', 'ENST0002|SNV-200-G-C')
]

class TestDecoyFasta(TestCaseIntegration, TestFastaWriterMixin):
    """ Test cases for decoyFasta """
    def generate_default_args(self):
        """ Generate the default args """
        args = argparse.Namespace()
        args.command = 'decoyFasta'
        args.output_path = self.work_dir/'test_output.fasta'
        args.decoy_string = 'DECOY_'
        args.decoy_string_position = 'prefix'
        args.method = 'reverse'
        args.shuffle_max_attempts = float('inf')
        args.non_shuffle_pattern = ''
        args.keep_peptide_nterm = 'true'
        args.keep_peptide_cterm = 'true'
        args.seed = None
        args.order = 'juxtaposed'
        args.quiet = False
        return args

    def test_decoy_fasta_case1(self):
        """ Test decoyFasta """
        self.write_test_fasta(TARGET_DB)
        args = self.generate_default_args()
        args.input_path = self.work_dir/'test_input.fasta'
        cli.decoy_fasta(args)

        received = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {args.input_path.name, args.output_path.name}
        self.assertEqual(received, expected)

        seqs = set(str(x.seq) for x in SeqIO.parse(args.output_path, 'fasta'))
        expected = {'MVPPILQWQREGGTPIK', 'ISSSCVCCIIQGGPPRR'}
        self.assertTrue(expected.issubset(seqs))

    def test_decoy_fasta_case2(self):
        """ Test decoyFasta """
        self.write_test_fasta(TARGET_DB)
        args = self.generate_default_args()
        args.method = 'shuffle'
        args.input_path = self.work_dir/'test_input.fasta'
        args.seed = 123123
        cli.decoy_fasta(args)

        received = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {args.input_path.name, args.output_path.name}
        self.assertEqual(received, expected)

        seqs = list(SeqIO.parse(args.output_path, 'fasta'))
        headers = {seq.description for seq in seqs}
        expected = {
            'ENST0001|SNV-100-A-T', 'DECOY_ENST0001|SNV-100-A-T',
            'ENST0002|SNV-200-G-C', 'DECOY_ENST0002|SNV-200-G-C'
        }
        self.assertEqual(headers, expected)

        seqs = {str(seq.seq) for seq in seqs}
        expected = {'MIPQVIWQEGGLPTRPK', 'IQCCGCIISRPSVSGPR'}
        self.assertTrue(expected.issubset(seqs))

    def test_decoy_fasta_case3(self):
        """ Test decoyFasta """
        self.write_test_fasta(TARGET_DB)
        args = self.generate_default_args()
        args.method = 'shuffle'
        args.input_path = self.work_dir/'test_input.fasta'
        args.seed = 123123
        args.non_shuffle_pattern = 'K,R'
        cli.decoy_fasta(args)

        received = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {args.input_path.name, args.output_path.name}
        self.assertEqual(received, expected)

        seqs = {str(seq.seq) for seq in SeqIO.parse(args.output_path, 'fasta')}
        expected = {'MIPWVIQREGGLPTQPK', 'IRIVSQCCISGPPCGSR'}
        self.assertTrue(expected.issubset(seqs))

    def test_decoy_fasta_cli(self):
        """ Test decoyFasta CLI """
        cmd = f"""
        {sys.executable} -m moPepGen.cli decoyFasta \\
            -i {self.data_dir}/peptides/variant.fasta \\
            -o {self.work_dir}/test_output.fasta
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise
