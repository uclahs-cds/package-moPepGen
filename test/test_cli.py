""" Test the command line interface """
import unittest
import shutil
from pathlib import Path
import argparse
from Bio import SeqIO
from moPepGen import cli, seqvar


BASE_DIR = Path(__file__).parent.absolute()
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
        shutil.rmtree(WORK_DIR, ignore_errors=True)

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
        args.verbose = False
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'star_fusion.tvf'}
        self.assertEqual(files, expected)

    def test_parse_fusion_catcher(self):
        """ Test parseFusionCatcher """
        args = argparse.Namespace()
        args.fusion = DATA_DIR/'fusion/fusion_catcher.txt'
        args.index_dir = DATA_DIR/'index'
        args.output_prefix = str(WORK_DIR/'fusion_catcher')
        args.verbose = True
        cli.parse_fusion_catcher(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'fusion_catcher.tvf'}
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
        args.min_length = 7
        args.max_length = 25
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
        args.min_length = 7
        args.max_length = 25
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
        args.min_length = 7
        args.max_length = 25
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
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_varaint_peptide_case5(self):
        """ A test case with reported in issue #25, with 3 indel. """
        args = argparse.Namespace()
        args.input_variant = [
            str(DATA_DIR/'vep/ENST00000308182.9_CPCG0100_indel.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = WORK_DIR/'vep_moPepGen.fasta'
        args.index_dir = DATA_DIR/'downsampled_index/ENST00000308182.9'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in WORK_DIR.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)
        peptides = list(SeqIO.parse(WORK_DIR/'vep_moPepGen.fasta', 'fasta'))

        seqs = {str(seq.seq) for seq in peptides}
        expected = {
            'EVVKGPGAPPPGK',
            'EVVKGPGAPPPGKR',
            'EVVKGPVLMPHFPPGK',
            'EVVKGPVLMPHFPPGKR',
            'EVVKGPVPHLER',
            'EVVQGSSAPAASS',
            'EVVQGSSAPAASSPTWK',
            'EVVQGSSAPAASSPTWKEVVK',
            'EVVQGSSAPAASVLLMPHFPPGK',
            'EVVQGSSAPAASVLLMPHFPPGKR',
            'EVVQGSSAPAASVLLMPHFPPGKRL',
            'EVVQGSSAPAASVLPHLER',
            'GCEGPWCSCCLIAHPER',
            'GPGAPPPGK',
            'GPGAPPPGKR',
            'GPGAPPPGKRL',
            'GPVLMPHFPPGK',
            'GPVLMPHFPPGKR',
            'GPVLMPHFPPGKRL',
            'GPVPHLER',
            'RGCEGPWCSCCLIAHPER'
        }
        self.assertEqual(seqs, expected)

    def test_parse_rmats_se_case_1(self):
        """ rMATS skipped exon when the retained version is annotated. This
        should results an deletion. """
        args = argparse.Namespace()
        args.skipped_exon = Path('test/files/alternative_splicing/rmats_se_case_1.txt')
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.location.start, 323)
        self.assertTrue(record.id, 'SE_324')
        self.assertTrue(int(record.attrs['START']), 323)
        self.assertTrue(int(record.attrs['END']), 405)
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_se_case_2(self):
        """ rMATS skipped exon when the skipped version is annotated. This
        should result an insertion. """
        args = argparse.Namespace()
        args.skipped_exon = Path('test/files/alternative_splicing/rmats_se_case_2.txt')
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.location.start, 870)
        self.assertTrue(record.id, 'SE_870')
        self.assertTrue(int(record.attrs['START']), 870)
        self.assertTrue(int(record.attrs['END']), 1097)
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_a5ss_case_1(self):
        """ rMATS A5SS when the longer version is annotated. This should
        results a deletion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = Path('test/files/alternative_splicing/rmats_a5ss_case_1.txt')
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_a5ss_case_2(self):
        """ rMATS A5SS when the shorter version is annotated. This should
        results an insertion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = Path('test/files/alternative_splicing/rmats_a5ss_case_2.txt')
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_a3ss_case_1(self):
        """ rMATS A3SS when the longer version is annotated. This should
        results a deletion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = Path('test/files/alternative_splicing/rmats_a3ss_case_1.txt')
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_a3ss_case_2(self):
        """ rMATS A3SS when the shorter version is annotated. This should
        results an Insertion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = Path('test/files/alternative_splicing/rmats_a3ss_case_2.txt')
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_mxe_case_1(self):
        """ rMATS MXE when one exon is annotated. This should results a
        substitution. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = Path('test/files/alternative_splicing/rmats_mxe_case_1.txt')
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Substitution')

    def test_parse_rmats_mxe_case_2(self):
        """ rMATS MXE when both exons are annotated. This should results two
        deletions. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = Path('test/files/alternative_splicing/rmats_mxe_case_2.txt')
        args.retained_intron = None
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))
        self.assertTrue(len(records), 2)
        for record in records:
            self.assertEqual(record.type, 'Deletion')

    def test_parse_rmats_ri(self):
        """ rMATS RI. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = Path('test/files/alternative_splicing/rmats_ri_case_1.txt')
        args.index_dir = Path('test/files/index')
        args.output_prefix = str(WORK_DIR/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))
        self.assertTrue(len(records), 2)
        for record in records:
            self.assertEqual(record.type, 'Insertion')


if __name__ == '__main__':
    unittest.main()
