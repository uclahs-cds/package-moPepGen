""" Integration test of splitFasta """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestSummarizeFasta(TestCaseIntegration):
    """ Test cases for summarizeFasta """
    def create_base_args(self) -> argparse.Namespace:
        """  Create base args """
        args = argparse.Namespace()
        args.command = 'summarizeFasta'
        args.quiet = False
        args.index_dir = None
        args.order_source = None
        args.output_path = self.work_dir/'output.txt'
        return args

    def test_summarize_fasta_case1(self):
        """ summarize fasta case1 """
        args = self.create_base_args()
        args.gvf = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'vep/vep_gINDEL.gvf',
            self.data_dir/'reditools/reditools.gvf',
            self.data_dir/'fusion/star_fusion.gvf',
            self.data_dir/'circRNA/circ_rna.gvf'
        ]
        args.variant_peptides = self.data_dir/'peptides/variant.fasta'
        args.annotation_gtf = self.data_dir/"annotation.gtf"
        cli.summarize_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {args.output_path.name}
        self.assertEqual(files, expected)

    def test_summarize_fasta_cli(self):
        """ Test summarizeFasta cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli summarizeFasta \\
            --gvf \\
                {self.data_dir}/vep/vep_gSNP.gvf \\
                {self.data_dir}/vep/vep_gINDEL.gvf \\
                {self.data_dir}/reditools/reditools.gvf \\
                {self.data_dir}/fusion/star_fusion.gvf \\
                {self.data_dir}/circRNA/circ_rna.gvf \\
            --variant-peptides {self.data_dir}/peptides/variant.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            -o {self.work_dir}/test
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise