""" Test the CLI for callNoncoding """
import argparse
import os
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli


def create_base_args() -> argparse.Namespace:
    """ Create a base args """
    args = argparse.Namespace()
    args.command = 'callNoncodingPeptide'
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.reference_source = None
    args.output_path = None
    args.output_orf = None
    args.inclusion_biotypes = None
    args.exclusion_biotypes = None
    args.min_tx_length = 21
    args.orf_assignment = 'max'
    args.w2f_reassignment = False
    args.cleavage_rule = 'trypsin'
    args.cleavage_exception = 'trypsin_exception'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.quiet = True
    return args

class TestCallNoncodingPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callNoncoding """

    def test_call_noncoding_cli(self):
        """ test callNoncoding cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli callNoncoding \\
            -o {self.work_dir}/circ.fasta \\
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

    def test_call_noncoding_peptides_case1(self):
        """ test call noncoding peptides """
        args = create_base_args()
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        cli.call_noncoding_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding_peptide.fasta', 'noncoding_orf.fasta'}
        self.assertEqual(files, expected)
        with open(self.work_dir/'noncoding_peptide.fasta', 'rt') as handle:
            peptides = list(SeqIO.parse(handle, 'fasta'))
            ids = [p.id for p in peptides]
            self.assertEqual(len(ids),len(set(ids)))
            self.assertTrue(peptides[0].id.split('|')[0] != 'None')

    def test_call_noncoding_peptides_case2(self):
        """ test call noncoding peptides when no ORF is found """
        args = create_base_args()
        ref_dir = self.data_dir/'downsampled_reference/ENST00000644482.1'
        args.genome_fasta = ref_dir/'genome.fasta'
        args.annotation_gtf = ref_dir/'annotation.gtf'
        args.proteome_fasta = ref_dir/'proteome.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        cli.call_noncoding_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding_peptide.fasta', 'noncoding_orf.fasta'}
        self.assertEqual(files, expected)
        size = os.stat(self.work_dir/'noncoding_peptide.fasta').st_size
        self.assertEqual(size, 0)
        size = os.stat(self.work_dir/'noncoding_orf.fasta').st_size
        self.assertEqual(size, 0)

    def test_call_noncoding_peptides_w2f(self):
        """ With w2f reassignment """
        args = create_base_args()
        args.w2f_reassignment = True
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        cli.call_noncoding_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding_peptide.fasta', 'noncoding_orf.fasta'}
        self.assertEqual(files, expected)
        with open(self.work_dir/'noncoding_peptide.fasta', 'rt') as handle:
            peptides = list(SeqIO.parse(handle, 'fasta'))
            ids = [p.id for p in peptides]
            self.assertEqual(len(ids),len(set(ids)))
            self.assertTrue(peptides[0].id.split('|')[0] != 'None')
            self.assertTrue(any('W2F' in x for x in ids))
