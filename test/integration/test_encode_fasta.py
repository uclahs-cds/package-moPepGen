""" Integration test for encodeFasta """
import argparse
from pathlib import Path
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestEncodeFasta(TestCaseIntegration):
    """ Test cases for encodeFasta """
    def test_encode_fasta(self):
        """ Test encodeFasta """
        args = argparse.Namespace()
        args.input_path = Path('test/files/peptides/variant.fasta')
        args.output_path = self.work_dir/'variant_encode.fasta'
        cli.encode_fasta(args)

        received = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'variant_encode.fasta', 'variant_encode.fasta.dict'}
        self.assertEqual(received, expected)

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
