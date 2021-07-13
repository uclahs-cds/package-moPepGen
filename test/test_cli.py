""" Test the command line interface """
import unittest
import shutil
import pathlib
import argparse
from moPepGen import cli


BASE_DIR = pathlib.Path(__file__).parent.absolute()
WORK_DIR = BASE_DIR / 'work_dir'
DATA_DIR = BASE_DIR / 'files'

class TestCli(unittest.TestCase):
    """ Test cases for the command line interface """
    def setUp(self):
        """ set up working directory """
        super().setUp()
        shutil.rmtree(WORK_DIR, ignore_errors=True)
        WORK_DIR.mkdir(parents=False, exist_ok=True)

    def tearDown(self):
        """ remove working files """
        super().tearDown()
        shutil.rmtree(WORK_DIR)

    def test_generate_index(self):
        """ Test genreate index """
        args = argparse.Namespace()
        args.genome_fasta = DATA_DIR / 'genome.fasta'
        args.annotation_gtf = DATA_DIR / 'annotation.gtf'
        args.proteome_fasta = DATA_DIR / 'translate.fasta'
        args.cleavage_rule = 'trypsin'
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.miscleavage = 2
        args.verbose = False
        args.output_dir = WORK_DIR / 'index'
        args.output_dir.mkdir(parents=False, exist_ok=True)
        cli.generate_index(args)
        files = {str(file.name) for file in args.output_dir.glob('*.pickle')}
        expected = {'genome.pickle', 'proteome.pickle', 'annotation.pickle',
            'canonical_peptides.pickle'}
        self.assertEqual(files, expected)

    def test_parse_vep(self):
        """ Test genreate index """
        args = argparse.Namespace()
        args.vep_txt = [
            DATA_DIR/'vep'/'vep_snp.txt',
            DATA_DIR/'vep'/'vep_indel.txt'
        ]
        # args.index_dir = None
        args.index_dir = DATA_DIR/'index'
        args.genome_fasta = DATA_DIR/'genome.fasta'
        args.annotation_gtf = DATA_DIR/'annotation.gtf'
        args.output_prefix = str(WORK_DIR/'vep')
        args.verbose = False
        cli.parse_vep(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep.tvf'}
        self.assertEqual(files, expected)

    def test_parse_star_fusion(self):
        """ Test parseSTARFusion """
        args = argparse.Namespace()
        args.fusion = DATA_DIR/'fusion/star_fusion.txt'
        args.index_dir = DATA_DIR/'index'
        args.output_prefix = str(WORK_DIR/'star_fusion')
        args.verbose = True
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'star_fusion.tvf'}
        self.assertEqual(files, expected)


    def test_call_variant_peptide_case1(self):
        """ Test variant peptide calling """
        args = argparse.Namespace()
        args.input_variant = [str(DATA_DIR/'vep'/'vep.tvf')]
        args.circ_rna_bed = None
        args.output_fasta = WORK_DIR/'vep_moPepGen.fasta'
        args.index_dir = DATA_DIR/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case2(self):
        """ Test variant peptide calling with fusion """
        args = argparse.Namespace()
        args.input_variant = [
            str(DATA_DIR/'vep'/'vep.tvf'),
            str(DATA_DIR/'fusion'/'fusion.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = WORK_DIR/'vep_moPepGen.fasta'
        args.index_dir = DATA_DIR/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case3(self):
        """ Test variant peptide calling with fusion and circRNA """
        args = argparse.Namespace()
        args.input_variant = [
            str(DATA_DIR/'vep'/'vep.tvf'),
            str(DATA_DIR/'fusion'/'fusion.tvf')
        ]
        args.circ_rna_bed = str(DATA_DIR/'circRNA'/'circ_rna.tsv')
        args.output_fasta = WORK_DIR/'vep_moPepGen.fasta'
        args.index_dir = DATA_DIR/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case4(self):
        """ Test variant peptide calling with alternative splicing """
        args = argparse.Namespace()
        args.input_variant = [
            str(DATA_DIR/'vep/vep.tvf'),
            str(DATA_DIR/'alternative_splicing/alternative_splicing.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = WORK_DIR/'vep_moPepGen.fasta'
        args.index_dir = DATA_DIR/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

if __name__ == '__main__':
    unittest.main()
