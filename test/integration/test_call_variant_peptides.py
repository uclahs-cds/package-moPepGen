""" Test the command line interface """
import argparse
from pathlib import Path
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli


class TestCallVariantPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callPeptides """

    def default_test_case(self, tvf:Path, reference:Path, expect:Path):
        """ Wrapper function to test actual cases.

        Args:
            tvf (Path): Path to the vep file.
            index (Path): Path to the index dir. This can be generated by the
                test/downsample_reference.py script.
            expect (Path): Path to the file of expected variant peptide
                sequence from the tvf file. This can be generated by the
                test/call_variant_peptide_brute_force.py script.
        """
        args = argparse.Namespace()
        args.input_variant = [str(tvf)]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.genome_fasta = reference/'genome.fasta'
        args.annotation_gtf = reference/'annotation.gtf'
        args.proteome_fasta = reference/'proteome.fasta'
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
        with open(expect, 'rt') as handle:
            expected = {line.strip() for line in handle}
        self.assertEqual(seqs, expected)

    def test_call_variant_peptide_case1(self):
        """ Test variant peptide calling """
        args = argparse.Namespace()
        args.input_variant = [str(self.data_dir/'vep'/'vep.tvf')]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
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
            str(self.data_dir/'vep'/'vep.gvf'),
            str(self.data_dir/'fusion'/'fusion.gvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
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
            str(self.data_dir/'vep'/'vep.gvf'),
            str(self.data_dir/'fusion'/'fusion.tvf')
        ]
        args.circ_rna_bed = str(self.data_dir/'circRNA'/'circ_rna.tsv')
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
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
            str(self.data_dir/'alternative_splicing/alternative_splicing.gvf')
        ]
        args.circ_rna_bed = None
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
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
        """ A test case reported in issue #25, with 3 indel. """
        tvf = self.data_dir \
            /'vep/CPCG0100_gencode_aa_indel_ENST00000308182.9.gvf'
        expect = self.data_dir \
            /'vep/CPCG0100_gencode_aa_indel_ENST00000308182.9_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000308182.9'
        self.default_test_case(tvf, reference, expect)

    def test_call_varaint_peptide_case6(self):
        """ A test case reported in issue #33, with 3 indel (insertion).
        """
        tvf = self.data_dir \
            /'vep/CPCG0102_gencode_aa_indel_ENST00000542218.1.gvf'
        expect = self.data_dir \
            /'vep/CPCG0102_gencode_aa_indel_ENST00000542218.1_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000542218.1'
        self.default_test_case(tvf, reference, expect)

    def test_call_varaint_peptide_case7(self):
        """ A test case reported in issue #33, with 3 indel (insertion).
        """
        tvf = self.data_dir \
            /'vep/CPCG0103_gencode_aa_indel_ENST00000314675.11.gvf'
        expect = self.data_dir \
            /'vep/CPCG0103_gencode_aa_indel_ENST00000314675.11_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(tvf, reference, expect)

    def test_call_varaint_peptide_case8(self):
        """ A test case reported in PR #36.
        """
        tvf = self.data_dir \
            /'vep/CPCG0184_gencode_aa_indel_ENST00000314675.11.gvf'
        expect = self.data_dir \
            /'vep/CPCG0184_gencode_aa_indel_ENST00000314675.11_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(tvf, reference, expect)

    def test_call_varaint_peptide_case9(self):
        """ A test case reported in PR #36.
        """
        tvf = self.data_dir \
            /'vep/CPCG0102_gencode_aa_indel_ENST00000360004.5.gvf'
        expect = self.data_dir \
            /'vep/CPCG0102_gencode_aa_indel_ENST00000360004.5_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000360004.5'
        self.default_test_case(tvf, reference, expect)

    def test_call_varaint_peptide_case10(self):
        """ A test case reported in PR #36.
        """
        tvf = self.data_dir \
            /'vep/CPCG0361_gencode_aa_indel_ENST00000390283.2.gvf'
        expect = self.data_dir \
            /'vep/CPCG0361_gencode_aa_indel_ENST00000390283.2_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000390283.2'
        self.default_test_case(tvf, reference, expect)

    def test_call_varaint_peptide_case11(self):
        """ A test case reported in PR #36. Testing the variant in the end
        of the sequence is included.
        """
        tvf = self.data_dir \
            /'vep/CPCG0100_gencode_aa_indel_ENST00000515757.5.tvf'
        expect = self.data_dir \
            /'vep/CPCG0100_gencode_aa_indel_ENST00000515757.5_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000515757.5'
        self.default_test_case(tvf, reference, expect)
