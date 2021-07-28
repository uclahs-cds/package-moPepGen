""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli


class TestCallVariantPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callPeptides """

    def test_call_variant_peptide_case1(self):
        """ Test variant peptide calling """
        args = argparse.Namespace()
        args.input_variant = [str(self.data_dir/'vep'/'vep.tvf')]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case2(self):
        """ Test variant peptide calling with fusion """
        args = argparse.Namespace()
        args.input_variant = [
            str(self.data_dir/'vep'/'vep.tvf'),
            str(self.data_dir/'fusion'/'fusion.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case3(self):
        """ Test variant peptide calling with fusion and circRNA """
        args = argparse.Namespace()
        args.input_variant = [
            str(self.data_dir/'vep'/'vep.tvf'),
            str(self.data_dir/'fusion'/'fusion.tvf')
        ]
        args.circ_rna_bed = str(self.data_dir/'circRNA'/'circ_rna.tsv')
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case4(self):
        """ Test variant peptide calling with alternative splicing """
        args = argparse.Namespace()
        args.input_variant = [
            str(self.data_dir/'vep/vep.tvf'),
            str(self.data_dir/'alternative_splicing/alternative_splicing.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'index'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_varaint_peptide_case5(self):
        """ A test case with reported in issue #25, with 3 indel. """
        args = argparse.Namespace()
        args.input_variant = [
            str(self.data_dir/'vep/CPCG0100_gencode_aa_indel_ENST00000308182.9.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'downsampled_index/ENST00000308182.9'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)
        peptides = list(SeqIO.parse(self.work_dir/'vep_moPepGen.fasta', 'fasta'))

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

    def test_call_varaint_peptide_case6(self):
        """ A test case with reported in issue #33, with 3 indel (insertion).
        """
        args = argparse.Namespace()
        args.input_variant = [
            str(self.data_dir/'vep/CPCG0102_gencode_aa_indel_ENST00000542218.1.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'downsampled_index/ENST00000542218.1'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)
        peptides = list(SeqIO.parse(self.work_dir/'vep_moPepGen.fasta', 'fasta'))
        seqs = {str(seq.seq) for seq in peptides}
        expected = {
            'EHTCQVMK',
            'EHTCQVMKR',
            'EHTCQVMKRHPN',
            'ERPWPREHTCQVMK',
            'ERPWPREHTCQVMKR',
            'ERPWPRGAHLSGYEEAS',
            'ERPWPRGEHLSGYEEAS',
            'ERPWPRGEHTCQVMK',
            'ERPWPRGEHTCQVMKR',
            'ERPWPRGGAHLSGYEEAS',
            'ERPWPRGGTPVR',
            'ERPWPRGGTPVRL',
            'ERPWPRGTPVR',
            'ERPWPRGTPVRL',
            'GAHLSGYEEAS',
            'GEHLSGYEEAS',
            'GEHTCQVMK',
            'GEHTCQVMKR',
            'GEHTCQVMKRHPN',
            'GGAHLSGYEEAS',
            'GGTPVRL',
            'GQKERPWPREHTCQVMK',
            'GQKERPWPRGAHLSGYEEAS',
            'GQKERPWPRGEHLSGYEEAS',
            'GQKERPWPRGEHTCQVMK',
            'GQKERPWPRGGAHLSGYEEAS',
            'GQKERPWPRGGTPVR',
            'GQKERPWPRGTPVR'
        }
        self.assertEqual(seqs, expected)

    def test_call_varaint_peptide_case7(self):
        """ A test case with reported in issue #33, with 3 indel (insertion).
        """
        args = argparse.Namespace()
        args.input_variant = [
            str(self.data_dir/'vep/CPCG0103_gencode_aa_indel_ENST00000314675.11.tvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = self.data_dir/'downsampled_index/ENST00000314675.11'
        args.cleavage_rule = 'trypsin'
        args.miscleavage = '2'
        args.min_mw = '500.'
        args.min_length = 7
        args.max_length = 25
        args.verbose = False
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)
        peptides = list(SeqIO.parse(self.work_dir/'vep_moPepGen.fasta', 'fasta'))
        seqs = {str(seq.seq) for seq in peptides}
        expected = {
            'APAPSTR',
            'APAPSTRCSAR',
            'APAPSTRCSARLLGGPSR',
            'APKSSLKFSPGPCPGPGPGPSPSR',
            'APKSSLKFSPGPCPGPGPSPSR',
            'AWDQLGQGGGGADPPLGNSSPGWR',
            'AWLLGLCLLGSSFSQEGLSGSCEGR',
            'CSARLLGGPSR',
            'CSARLLGGPSRQAR',
            'FSPGPCPGPGPGPSPSR',
            'FSPGPCPGPGPGPSPSRPQSR',
            'FSPGPCPGPGPGPSPSRPQSRSR',
            'FSPGPCPGPGPSPSR',
            'FSPGPCPGPGPSPSRPQSR',
            'FSPGPCPGPGPSPSRPQSRSR',
            'GLGPAGPGR',
            'GLGPAGPGRGR',
            'GLGPAGPGRGRCRPAPW',
            'GRCRPAPW',
            'HPPPQPAAPPGSWVGHLAR',
            'HPPPQPAAPPGSWVGHLARHVGK',
            'HRPTNSWGHFNTAQSHFATRPMPPN',
            'KGLGPAGPGR',
            'KGLGPAGPGRGR',
            'LLGGPSR',
            'LLGGPSRQAR',
            'LLGGPSRQARR',
            'LPAPVPVLDPVPAPNK',
            'LPAPVPVLDPVPAPNKAPAPSTR',
            'LQSLSWTQSQPPIK',
            'RKGLGPAGPGR',
            'SRPQSLSWTQSQPPIK',
            'SRSRPQSLSWTQSQPPIK',
            'SRSSPCPGPSPSPQ',
            'SSLKFSPGPCPGPGPGPSPSR',
            'SSLKFSPGPCPGPGPGPSPSRPQSR',
            'SSLKFSPGPCPGPGPSPSR',
            'SSLKFSPGPCPGPGPSPSRPQSR',
            'SSLKFSPGPCPGPGPSPSRPQSRSR',
            'SSPCPGPSPSPQ'
        }
        self.assertEqual(seqs, expected)
