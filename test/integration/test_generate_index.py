""" Test the command line interface """
import argparse
import subprocess as sp
import os
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli, aa
from moPepGen.gtf import GenomicAnnotation, GenomicAnnotationOnDisk


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
        args.gtf_symlink = False
        args.reference_source = None
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
        files = {str(file.name) for file in args.output_dir.glob('*')}
        expected = {'genome.pkl', 'proteome.pkl',
            'annotation.gtf', 'annotation_gene.idx', 'annotation_tx.idx',
            'canonical_peptides.pkl', 'coding_transcripts.pkl'}
        self.assertEqual(files, expected)

class TestCaseGenomicAnnotationOnDisk(TestCaseIntegration):
    """ Test case for GenomicAnnotationOnDisk """
    def test_generate_index(self):
        """ Test generate index """
        proteome = aa.AminoAcidSeqDict()
        proteome.dump_fasta(self.data_dir/'translate.fasta')
        gtf_file = self.data_dir/'annotation.gtf'
        gtf_file2 = self.work_dir/'annotation.gtf'
        os.symlink(gtf_file, gtf_file2)
        cli.index_gtf(gtf_file2, proteome=proteome)
        expect = {
            'annotation.gtf', 'annotation_gene.idx', 'annotation_tx.idx'
        }
        received = {str(file.name) for file in self.work_dir.glob('*')}
        self.assertEqual(received, expect)


    def test_load_index(self):
        """ test load index """
        anno1 = GenomicAnnotation()
        anno1.dump_gtf(self.data_dir/'annotation.gtf')

        anno2 = GenomicAnnotationOnDisk()
        anno2.load_index(self.data_dir/'annotation.gtf')

        self.assertEqual(anno1.genes.keys(), anno2.genes.keys())
        self.assertEqual(anno1.transcripts.keys(), anno2.transcripts.keys())
