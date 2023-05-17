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
        args.cleavage_rule = 'trypsin'
        args.output_path = self.work_dir/'output.txt'
        args.output_image = None
        args.ignore_missing_source = False
        args.reference_source = None
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
        args.noncoding_peptides = self.data_dir/'peptides/noncoding.fasta'
        args.alt_translation_peptides = self.data_dir/'peptides/alt_translation.fasta'
        args.annotation_gtf = self.data_dir/"annotation.gtf"
        args.proteome_fasta = self.data_dir/"translate.fasta"
        cli.summarize_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {args.output_path.name}
        self.assertEqual(files, expected)

    def test_summarize_fasta_case2(self):
        """ summarize fasta case2 with sources missing """
        args = self.create_base_args()
        args.gvf = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'vep/vep_gINDEL.gvf',
            self.data_dir/'reditools/reditools.gvf',
            self.data_dir/'fusion/star_fusion.gvf',
            self.data_dir/'circRNA/circ_rna.gvf'
        ]
        args.variant_peptides = self.data_dir/'peptides/variant.fasta'
        args.noncoding_peptides = self.data_dir/'peptides/noncoding.fasta'
        args.alt_translation_peptides = None
        args.annotation_gtf = self.data_dir/"annotation.gtf"
        args.proteome_fasta = self.data_dir/"translate.fasta"
        args.order_source = 'gSNP,gINDEL,RNAEditingSite,Fusion,circRNA,rMATS,Noncoding'
        args.ignore_missing_source = True
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
            -p {self.data_dir}/translate.fasta \\
            -o {self.work_dir}/test.txt
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise
