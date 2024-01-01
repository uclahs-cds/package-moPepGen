""" Integration test of splitFasta """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestSplitDatabase(TestCaseIntegration):
    """ Test cases for splitFasta """
    def test_split_fasta_cli(self):
        """ test splitFasta cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli splitFasta \\
            --gvf \\
                {self.data_dir}/vep/vep_gSNP.gvf \\
                {self.data_dir}/vep/vep_gINDEL.gvf \\
                {self.data_dir}/reditools/reditools.gvf \\
                {self.data_dir}/fusion/star_fusion.gvf \\
                {self.data_dir}/circRNA/circ_rna.gvf \\
            --variant-peptides {self.data_dir}/peptides/variant.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta \\
            -o {self.work_dir}/test
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
        args.command = 'splitFasta'
        args.noncoding_peptides = None
        args.order_source = None
        args.group_source = None
        args.max_source_groups = 1
        args.additional_split = None
        args.output_prefix = self.work_dir/'test'
        args.reference_source = None
        args.index_dir = None
        args.quiet = True
        return args

    def test_split_fasta_case1(self):
        """ test splitFasta case1 """
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
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.split_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {
            'test_gINDEL.fasta','test_gSNP.fasta', 'test_Fusion.fasta',
            'test_RNAEditingSite.fasta', 'test_circRNA.fasta',
            'test_Remaining.fasta', 'test_circRNA.fasta', 'test_Noncoding.fasta',
            'test_CodonReassign.fasta', 'test_SECT.fasta'
        }
        self.assertEqual(files, expected)

    def test_split_fasta_case2(self):
        """ test splitFasta with sources being grouped """
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
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.group_source = ['coding:gSNP,gINDEL']
        cli.split_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {
            'test_coding.fasta', 'test_RNAEditingSite.fasta', 'test_Fusion.fasta',
            'test_circRNA.fasta', 'test_Remaining.fasta', 'test_circRNA.fasta',
            'test_Noncoding.fasta', 'test_CodonReassign.fasta', 'test_SECT.fasta'
        }
        self.assertEqual(files, expected)

    def test_split_fasta_case3(self):
        """ test splitFasta with groups and order. """
        args = self.create_base_args()
        args.gvf = [
            self.data_dir/'reditools/reditools.gvf',
            self.data_dir/'fusion/star_fusion.gvf',
            self.data_dir/'circRNA/circ_rna.gvf',
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'vep/vep_gINDEL.gvf'
        ]
        args.variant_peptides = self.data_dir/'peptides/variant.fasta'
        args.noncoding_peptides = self.data_dir/'peptides/noncoding.fasta'
        args.alt_translation_peptides = self.data_dir/'peptides/alt_translation.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.group_source = [
            'DNA:gSNP,gINDEL',
            'RNA:RNAEditingSite,Fusion,circRNA',
            'ALT:SECT,CodonReassign'
        ]
        args.order_source = 'ALT,DNA,RNA,Noncoding'
        cli.split_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {
            'test_DNA.fasta', 'test_RNA.fasta',
            'test_ALT.fasta', 'test_Noncoding.fasta',
            'test_Remaining.fasta'
        }
        self.assertEqual(files, expected)

    def test_split_fasta_case4(self):
        """ test splitFasta case 4 with additional split """
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
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.additional_split = ['Noncoding-gSNP']
        cli.split_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {
            'test_gINDEL.fasta','test_gSNP.fasta', 'test_Fusion.fasta',
            'test_RNAEditingSite.fasta', 'test_circRNA.fasta',
            'test_Remaining.fasta', 'test_circRNA.fasta', 'test_Noncoding.fasta',
            'test_CodonReassign.fasta', 'test_SECT.fasta'
        }
        self.assertEqual(files, expected)

    def test_split_fasta_case5(self):
        """ test splitFasta altTrans are grouped """
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
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.group_source = ['coding:gSNP,gINDEL', 'ALT:SECT,CodonReassign']
        cli.split_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {
            'test_coding.fasta', 'test_RNAEditingSite.fasta', 'test_Fusion.fasta',
            'test_circRNA.fasta', 'test_Remaining.fasta', 'test_circRNA.fasta',
            'test_Noncoding.fasta', 'test_ALT.fasta'
        }
        self.assertEqual(files, expected)

    def test_split_fasta_source_order_comb(self):
        """ test splitFasta with source order of combinations """
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
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.group_source = [
            'ALT:SECT,CodonReassign',
            'NotCirc:gSNP,gINDEL,sSNV,sINDEL,Fusion,altSplice,RNAEditingSite'
        ]
        args.order_source = ','.join([
            'NotCir',
            'ALT',
            'NotCirc-ALT',
            'Noncoding',
            'Noncoding-NotCirc',
            'Noncoding-ALT',
            'circRNA',
            'circRNA-ALT',
            'circRNA-NotCirc',
            'Noncoding-circRNA'
        ])
        args.max_source_groups = 2
        cli.split_fasta(args)

    def test_split_fasta_noncoding_and_alttrans(self):
        """ test splitFasta with only noncoding and alt trans peptides """
        args = self.create_base_args()
        args.gvf = [
        ]
        args.variant_peptides = None
        args.noncoding_peptides = self.data_dir/'peptides/noncoding.fasta'
        args.alt_translation_peptides = self.data_dir/'peptides/alt_translation.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.split_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {
            'test_Noncoding.fasta', 'test_SECT.fasta',
            'test_Remaining.fasta', 'test_CodonReassign.fasta'
        }
        self.assertEqual(files, expected)
