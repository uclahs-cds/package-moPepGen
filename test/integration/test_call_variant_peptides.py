""" Test the command line interface """
import argparse
from typing import List
from pathlib import Path
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli


def create_base_args() -> argparse.Namespace:
    """ Create base args """
    args = argparse.Namespace()
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.output_fasta = None
    args.max_variants_per_node = 5
    args.inclusion_biotypes = None
    args.exclusion_biotypes = None
    args.cleavage_rule = 'trypsin'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.quiet = True
    args.verbose_level = 1
    args.noncanonical_transcripts = False
    args.invalid_protein_as_noncoding = False
    args.threads = 1
    return args

class TestCallVariantPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callPeptides """

    def default_test_case(self, gvf:List[Path], reference:Path, expect:Path):
        """ Wrapper function to test actual cases.

        Args:
            tvf (Path): Path to the vep file.
            index (Path): Path to the index dir. This can be generated by the
                test/downsample_reference.py script.
            expect (Path): Path to the file of expected variant peptide
                sequence from the tvf file. This can be generated by the
                test/call_variant_peptide_brute_force.py script.
        """
        args = create_base_args()
        args.input_variant = gvf
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = reference/'genome.fasta'
        args.annotation_gtf = reference/'annotation.gtf'
        args.proteome_fasta = reference/'proteome.fasta'
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
        args = create_base_args()
        args.input_variant = [self.data_dir/'vep'/'vep_gSNP.gvf']
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case2(self):
        """ Test variant peptide calling with fusion """
        args = create_base_args()
        args.input_variant = [
            self.data_dir/'vep'/'vep_gSNP.gvf',
            self.data_dir/'fusion'/'fusion.gvf'
        ]
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case3(self):
        """ Test variant peptide calling with fusion and circRNA """
        args = create_base_args()
        args.input_variant = [
            self.data_dir/'vep'/'vep_gSNP.gvf',
            self.data_dir/'fusion'/'fusion.gvf',
            self.data_dir/'circRNA'/'circ_rna.gvf'
        ]
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_all_sources(self):
        """ Test variant peptide calling with fusion and circRNA,
        RNAEditing, gSNP and gINDEL """
        args = create_base_args()
        args.input_variant = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'vep/vep_gINDEL.gvf',
            self.data_dir/'fusion/fusion.gvf',
            self.data_dir/'circRNA/circ_rna.gvf',
            self.data_dir/'reditools/reditools.gvf'
        ]
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_case4(self):
        """ Test variant peptide calling with alternative splicing """
        args = create_base_args()
        args.input_variant = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'alternative_splicing/alternative_splicing.gvf'
        ]
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_varaint_peptide_case5(self):
        """ A test case reported in issue #25, with 3 indel. """
        gvf = [
            self.data_dir/'vep/CPCG0100_gencode_aa_indel_ENST00000308182.9.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0100_gencode_aa_indel_ENST00000308182.9_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000308182.9'
        self.default_test_case(gvf, reference, expect)

    def test_call_varaint_peptide_case6(self):
        """ A test case reported in issue #33, with 3 indel (insertion).
        """
        gvf = [
            self.data_dir/'vep/CPCG0102_gencode_aa_indel_ENST00000542218.1.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0102_gencode_aa_indel_ENST00000542218.1_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000542218.1'
        self.default_test_case(gvf, reference, expect)

    def test_call_varaint_peptide_case7(self):
        """ A test case reported in issue #33, with 3 indel (insertion).
        """
        gvf = [
            self.data_dir/'vep/CPCG0103_gencode_aa_indel_ENST00000314675.11.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0103_gencode_aa_indel_ENST00000314675.11_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expect)

    def test_call_varaint_peptide_case8(self):
        """ A test case reported in PR #36.
        """
        gvf = [
            self.data_dir/'vep/CPCG0184_gencode_aa_indel_ENST00000314675.11.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0184_gencode_aa_indel_ENST00000314675.11_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expect)

    def test_call_varaint_peptide_case9(self):
        """ A test case reported in PR #36.
        """
        gvf = [
            self.data_dir/'vep/CPCG0102_gencode_aa_indel_ENST00000360004.5.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0102_gencode_aa_indel_ENST00000360004.5_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000360004.5'
        self.default_test_case(gvf, reference, expect)

    def test_call_varaint_peptide_case10(self):
        """ A test case reported in PR #36.
        """
        gvf = [
            self.data_dir/'vep/CPCG0361_gencode_aa_indel_ENST00000390283.2.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0361_gencode_aa_indel_ENST00000390283.2_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000390283.2'
        self.default_test_case(gvf, reference, expect)

    def test_call_varaint_peptide_case11(self):
        """ A test case reported in PR #36. Testing the variant in the end
        of the sequence is included.
        """
        gvf = [
            self.data_dir/'vep/CPCG0100_gencode_aa_indel_ENST00000515757.5.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0100_gencode_aa_indel_ENST00000515757.5_expect.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000515757.5'
        self.default_test_case(gvf, reference, expect)

    def test_call_variant_peptide_case12(self):
        """ A test case that failed to find any variant peptide """
        gvf = [
            self.data_dir/'vep/CPCG0100_indel_ENST00000317799.10.gvf'
        ]
        expect = self.data_dir \
            /'vep/CPCG0100_indel_ENST00000317799.10_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000317799.10'
        self.default_test_case(gvf, reference, expect)

    def test_call_variant_peptide_case13(self):
        """ Reported in issue #165 """
        gvf = [
            self.data_dir/'vep/CPCG0339_indel_ENST00000360004.5.gvf'
        ]
        expected = self.data_dir \
            /'vep/CPCG0339_indel_ENST00000360004.5_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000360004.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case14(self):
        """ Reported in issue # """
        gvf = [
            self.data_dir/'vep/CPCG0233_indel_ENST00000314675.11.gvf'
        ]
        expected = self.data_dir \
            /'vep/CPCG0233_indel_ENST00000314675.11_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case15(self):
        """ Reported in issue #171 """
        gvf = [
            self.data_dir/'vep/CPCG0190_indel_ENST00000281589.4.gvf'
        ]
        expected = self.data_dir \
            /'vep/CPCG0190_indel_ENST00000281589.4_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000281589.4'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case16(self):
        """ Reported in issue #231 noncoding transcript """
        gvf = [
            self.data_dir/'vep/CPCG0100_ENST00000446393.2.gvf'
        ]
        expected = self.data_dir \
            /'vep/CPCG0100_ENST00000446393.2_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000446393.2'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case17(self):
        """ Another noncoding transcript """
        gvf = [
            self.data_dir/'vep/CPCG0102_ENST00000485710.5.gvf'
        ]
        expected = self.data_dir \
            /'vep/CPCG0100_ENST00000485710.5_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000485710.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case18(self):
        """ snv + indel reported in #302. This test case has a region of its
        peptide graph that looks like below when creating the cleavage graph.
        The node 'MSKNIH' has two incoming nodes and two outgoing nodes. One
        of the incoming and one of the outgoing node is bridge to/from other
        reading frame. When the cursor is at node 'IIFFFWYLTLALNGQLFFFL', merge
        backwards will be applied and mergeing it with MSGNIH. The incoming 'L'
        of MSGNIH will also merge with it. The returned downstream nodes must
        be the merged 'IIFFFWYLTLALNGQLFFFLMSGNIH', rather than it's actual
        outgoing F and L, otherwise the downstream MNCRCLRNK*YSRARWLMLVI...IFI
        may remain unprocessed to cause * in the amino acid sequence, which
        causes the issue.

                                         /
                                        F
                                       /
        1:  IIFFFWYLTLALNGQLFFFL-MSGNIH-L-MNCRCLRNK*YSRARWLMLVI...IFI
                                /        /
                               L   *VEIFIL
                               |  / 219
        2:  -*FSFFGT*HWL*TDNFFFF-*-*VEIFIY-

        """
        gvf = [
            self.data_dir/'vep/CPCG0100_ENST00000623909.1/gindel.gvf',
            self.data_dir/'vep/CPCG0100_ENST00000623909.1/gsnp.gvf'
        ]
        expected = self.data_dir/'vep/CPCG0100_ENST00000623909.1_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000623909.1/'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case19(self):
        """ Fusion + variants reported in #299. Noted that the expected
        is generated by calling callVariant directly, because the bruteForce
        does not support fusion at this moment. Should update it later. """
        gvf = [
            self.data_dir/'comb/CPCG0100_ENST00000265138.4/arriba.gvf',
            self.data_dir/'comb/CPCG0100_ENST00000265138.4/gindel.gvf',
            self.data_dir/'comb/CPCG0100_ENST00000265138.4/gsnp.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0100_ENST00000265138.4_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptides_case20(self):
        """ Issue reported in #313. In this test case, some nodes are missed
        when creating the cleavage graph, and they remained uncleaved,
        resulting * in the peptide sequence. The solution is just to use
        `merge_forward` when the only outgoing node is a bridge, instead of
        `merge_join`. However, 5 peptides are sometimes not called. The graph
        is too complicated to figure out where is the cause. Going to leave it
        here for now.
        """
        gvf = [
            self.data_dir/'comb/CPCG0235_ENST00000480694.2/gsnp.gvf',
            self.data_dir/'comb/CPCG0235_ENST00000480694.2/gindel.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0235_ENST00000480694.2_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000480694.2'
        args = create_base_args()
        args.input_variant = gvf
        args.output_fasta = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = reference/'genome.fasta'
        args.annotation_gtf = reference/'annotation.gtf'
        args.proteome_fasta = reference/'proteome.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        self.assertEqual(files, {'vep_moPepGen.fasta'})
        peptides = list(SeqIO.parse(self.work_dir/'vep_moPepGen.fasta', 'fasta'))
        seqs = {str(seq.seq) for seq in peptides}
        with open(expected, 'rt') as handle:
            expect_seqs = {line.strip() for line in handle}
        self.assertTrue(seqs.issuperset(expect_seqs))

    def test_call_variant_peptide_case21(self):
        """ Alternative splicing only, reported in #333.
        NOTE: The expected sequence is created by running callVariant because
        bruteForce does not support for alternative splicing at this point.
        This is just to ensure callVariant runs on this situation.
        """
        gvf = [
            self.data_dir/'alternative_splicing/CPCG0486_ENST00000481806.1.gvf'
        ]
        expected = self.data_dir/'alternative_splicing/CPCG0486_ENST00000481806.1_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000481806.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case22(self):
        """ Fusion + gSNP, reported in #320.
        NOTE: The expected sequence is created by running callVariant
        """
        gvf = [
            self.data_dir/'comb/CPCG0464_ENST00000370143.5_ENST00000370165.7/arriba.gvf',
            self.data_dir/'comb/CPCG0464_ENST00000370143.5_ENST00000370165.7/gsnp.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0464_ENST00000370143.5_ENST00000370165.7_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000370143.5_ENST00000370165.7'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case23(self):
        """ SNP only variant from fuzz test. This test case is used to fix two
        issues:
        1. When calculating the orf stop site, it did not consider when it's
        on the negative strand
        2. The raised by stop codon is not included in the sequence when
        the sequence is truncated at the orf stop site. """
        gvf = [
            self.data_dir/'fuzz/01/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/01/brute_force.fasta'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case24(self):
        """ INDEL only variant from fuzz test. This test case ensures that
        when frameshifting mutations are incorporated to the TVG graph, not only
        the target reading frame from canonical reading frame is activated,
        but also from the other active reading frames. """
        gvf = [
            self.data_dir/'fuzz/02/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/02/brute_force.fasta'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case25(self):
        """ INDEL only variant from fuzz test. This test case ensures that
        when finding nodes to merge during creating the cleavage graph, all
        inbridge nodes are also included. """
        gvf = [
            self.data_dir/'fuzz/03/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/03/brute_force.fasta'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case26(self):
        """ Test case reported in #357 with two alternative splicing events. """
        gvf = [
            self.data_dir/'comb/CPCG0266_ENST00000381461.6/rMATs.gvf',
            self.data_dir/'comb/CPCG0266_ENST00000381461.6/gsnp.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0266_ENST00000381461.6_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000381461.6'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case27(self):
        """ Test case reported in #360 with one alternative splicing and an
        indel. """
        gvf = [
            self.data_dir/'comb/CPCG0235_ENST00000525687.5/rMATs.gvf',
            self.data_dir/'comb/CPCG0235_ENST00000525687.5/gindel.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0235_ENST00000525687.5_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000525687.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case28(self):
        """ Test case reported also in #360 with combination of one alternative
        splicing and indels & snps. The issue of the case is that when
        collapsing end nodes, the one with lower subgraph level should be kept.
        """
        gvf = [
            self.data_dir/'comb/CPCG0235_ENST00000590400.1/rMATs.gvf',
            self.data_dir/'comb/CPCG0235_ENST00000590400.1/gsnp.gvf',
            self.data_dir/'comb/CPCG0235_ENST00000590400.1/gindel.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0235_ENST00000590400.1_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000590400.1'
        self.default_test_case(gvf, reference, expected)
