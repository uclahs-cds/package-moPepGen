""" Integration test for encodeFasta """
import argparse
from pathlib import Path
import subprocess as sp
import sys
from test.integration import TestCaseIntegration, TestFastaWriterMixin
from Bio import SeqIO
from moPepGen import cli


class TestEncodeFasta(TestCaseIntegration, TestFastaWriterMixin):
    """ Test cases for encodeFasta """
    def generate_default_args(self):
        """ Generate default args """
        args = argparse.Namespace()
        args.command = 'encodeFasta'
        args.input_path = self.work_dir/'test_input.fasta'
        args.output_path = self.work_dir/'variant_encode.fasta'
        args.decoy_string = 'DECOY_'
        args.decoy_string_position = 'prefix'
        args.quiet = True
        return args

    def test_encode_fasta(self):
        """ Test encodeFasta """
        args = self.generate_default_args()
        args.input_path = Path('test/files/peptides/variant.fasta')
        cli.encode_fasta(args)

        received = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'variant_encode.fasta', 'variant_encode.fasta.dict'}
        self.assertEqual(received, expected)

    def test_encode_fasta_decoyed(self):
        """ Test encodeFasta with decoyed database """
        data = [
            ('MIPTGGERQWQLIPPVK', 'ENST0001|SNV-100-A-T'),
            ('KVPPILQWQREGGTPIM', 'DECOY_ENST0001|SNV-100-A-T'),
            ('IRPPGGQIICCVCSSSR', 'ENST0002|SNV-200-G-C'),
            ('RSSSCVCCIIQGGPPRI', 'DECOY_ENST0002|SNV-200-G-C')
        ]
        args = self.generate_default_args()
        self.write_test_fasta(data)
        cli.encode_fasta(args)
        received = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'variant_encode.fasta', 'variant_encode.fasta.dict', 'test_input.fasta'}
        self.assertTrue(received, expected)

        with open(self.work_dir/'test_input.fasta', 'rt') as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
        self.assertEqual(sum(x.description.startswith('DECOY_') for x in seqs), 2)
        self.assertEqual(len({x.description.replace('DECOY_', '') for x in seqs}), 2)

    def test_ecnode_fasta_cli(self):
        """ Test encodeFasta CLI """
        cmd = f"""
        {sys.executable} -m moPepGen.cli encodeFasta \\
            -i {self.data_dir}/peptides/variant.fasta \\
            -o {self.work_dir}/variant_encode.fasta
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise
