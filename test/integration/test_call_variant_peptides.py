""" Test the command line interface """
import argparse
from typing import List
from pathlib import Path
import sys
import subprocess as sp
from test.integration import TestCaseIntegration
from test.unit import create_variant
from Bio import SeqIO
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen import cli, seqvar
from moPepGen.aa import VariantPeptideIdentifier as vpi


def create_base_args() -> argparse.Namespace:
    """ Create base args """
    args = argparse.Namespace()
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.reference_source = None
    args.output_path = None
    args.max_variants_per_node = 7
    args.additional_variants_per_misc = 2
    args.min_nodes_to_collapse = 30
    args.naa_to_collapse = 5
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

    def default_test_case(self, gvf:List[Path], reference:Path, expect:Path=None):
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
        args.input_path = gvf
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = reference/'genome.fasta'
        args.annotation_gtf = reference/'annotation.gtf'
        args.proteome_fasta = reference/'proteome.fasta'
        args.reference_source = None
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)
        if not expect:
            return
        peptides = list(SeqIO.parse(self.work_dir/'vep_moPepGen.fasta', 'fasta'))
        seqs = {str(seq.seq) for seq in peptides}
        with open(expect, 'rt') as handle:
            expected = {line.strip() for line in handle}
        self.assertEqual(seqs, expected)

    def test_call_variant_cli(self):
        """ test callVariant cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli callVariant \\
            -i \\
                {self.data_dir}/vep/vep_gSNP.gvf \\
                    {self.data_dir}/vep/vep_gINDEL.gvf \\
                {self.data_dir}/fusion/fusion.gvf \\
                {self.data_dir}/circRNA/circ_rna.gvf \\
                {self.data_dir}/reditools/reditools.gvf \\
            -o {self.work_dir}/test.fasta \\
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

    def test_call_variant_peptide_case1(self):
        """ Test variant peptide calling """
        args = create_base_args()
        args.input_path = [self.data_dir/'vep'/'vep_gSNP.gvf']
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
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
        args.input_path = [
            self.data_dir/'vep'/'vep_gSNP.gvf',
            self.data_dir/'fusion'/'fusion.gvf'
        ]
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
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
        args.input_path = [
            self.data_dir/'vep'/'vep_gSNP.gvf',
            self.data_dir/'fusion'/'fusion.gvf',
            self.data_dir/'circRNA'/'circ_rna.gvf'
        ]
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
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
        args.input_path = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'vep/vep_gINDEL.gvf',
            self.data_dir/'fusion/fusion.gvf',
            self.data_dir/'circRNA/circ_rna.gvf',
            self.data_dir/'reditools/reditools.gvf'
        ]
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)

    def test_call_variant_peptide_fusion_only(self):
        """ Test variant peptide calling with fusion only. """
        args = create_base_args()
        args.input_path = [
            self.data_dir/'fusion'/'fusion.gvf'
        ]
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_moPepGen.fasta'}
        self.assertEqual(files, expected)
        seqs = list(SeqIO.parse(args.output_path, 'fasta'))
        self.assertTrue(len(seqs) > 0)

    def test_call_variant_peptide_fusion_donor_intronic(self):
        """ Test case that the donor breakpoint is intronic """
        attrs={
            'ACCEPTER_GENE_ID': 'ENSG00000244486.9',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST00000622235.5',
            'ACCEPTER_POSITION': '50',
            'TRANSCRIPT_ID': 'ENST00000642151.1'
        }
        variant = create_variant(
            start=94, end=95, ref='A', alt='<FUSION>', _type='Fusion',
            _id='ENST00000642151.1-94:ENST00000622235.5:50', attrs=attrs,
            seqname='ENSG00000099949.21'
        )
        metadata = GVFMetadata(
            parser='parseSTARFusion', source='Fusion', chrom='Gene ID'
        )
        seqvar.io.write([variant], self.work_dir/'fusion.gvf', metadata)
        args = create_base_args()
        args.input_path = [
            self.work_dir/'fusion.gvf'
        ]
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        cli.call_variant_peptide(args)
        seqs = list(SeqIO.parse(args.output_path, 'fasta'))
        self.assertTrue(len(seqs) > 0)

    def test_call_variant_peptide_case4(self):
        """ Test variant peptide calling with alternative splicing """
        args = create_base_args()
        args.input_path = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'alternative_splicing/alternative_splicing.gvf'
        ]
        args.output_path = self.work_dir/'vep_moPepGen.fasta'
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
        `merge_join`.
        """
        gvf = [
            self.data_dir/'comb/CPCG0235_ENST00000480694.2/gsnp.gvf',
            self.data_dir/'comb/CPCG0235_ENST00000480694.2/gindel.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0235_ENST00000480694.2_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000480694.2'
        self.default_test_case(gvf, reference, expected)

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

    def test_call_variant_peptide_case29(self):
        """ Test case reported in #364 with combination of two alternative
        splicing and indels & snps. The issue of the case is that when
        aligning nodes in TVG, it failed to limit inside the same subgraph.
        """
        gvf = [
            self.data_dir/'comb/CPCG0333_ENST00000452737.5/rMATs.gvf',
            self.data_dir/'comb/CPCG0333_ENST00000452737.5/gsnp.gvf',
            self.data_dir/'comb/CPCG0333_ENST00000452737.5/gindel.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0333_ENST00000452737.5_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case30(self):
        """ Test case reported in #364 with combination of two alternative
        splicing and indels & snps. The issue of the case is that when
        aligning nodes in TVG, it failed to limit inside the same subgraph.
        """
        gvf = [
            self.data_dir/'comb/CPCG0462_ENST00000483923.5/rMATs.gvf',
            self.data_dir/'comb/CPCG0462_ENST00000483923.5/gsnp.gvf'
        ]
        expected = self.data_dir/'comb/CPCG0462_ENST00000483923.5_expected.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000483923.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case31(self):
        """ Test case reported in #490 that variant label of the FASTA header
        was not written correctly when the entire variatn peptide is called
        from an insertion into a circRNA """
        gvf = [
            self.data_dir/'comb/CPCG0435_ENST00000374864.9/circ2.gvf',
            self.data_dir/'comb/CPCG0435_ENST00000374864.9/gsnp.gvf',
            self.data_dir/'comb/CPCG0435_ENST00000374864.9/gindel.gvf'
        ]
        reference = self.data_dir/'downsampled_reference/ENST00000374864.9'
        self.default_test_case(gvf, reference, None)
        has_incorrect_fasta_header = False
        for peptide in SeqIO.parse(self.work_dir/'vep_moPepGen.fasta', 'fasta'):
            labels = vpi.parse_variant_peptide_id(peptide.description)
            for label in labels:
                if not isinstance(label, vpi.CircRNAVariantPeptideIdentifier):
                    continue
                for var_id in label.variant_ids:
                    if var_id.startswith('ENS'):
                        has_incorrect_fasta_header = True
                        break
        self.assertFalse(has_incorrect_fasta_header)

    def test_call_variant_peptide_case32(self):
        """ Noncoding TX with start gain mutation from fuzz test. """
        gvf = [
            self.data_dir/'fuzz/04/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/04/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case33(self):
        """ Coding TX with deletion that start at forth nucleotide from CDS
        start site. """
        gvf = [
            self.data_dir/'fuzz/05/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/05/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case34(self):
        """ Noncoding TX with deletion that start at forth nucleotide from CDS
        start site. """
        gvf = [
            self.data_dir/'fuzz/06/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/06/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case35(self):
        """ Noncoding TX with frameshifting mutations. #508 #509 """
        gvf = [
            self.data_dir/'fuzz/07/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/07/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case36(self):
        """ Coding TX with in-frame deletion. #515 """
        gvf = [
            self.data_dir/'fuzz/08/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/08/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case37(self):
        """ Noncoding TX with stop lost mutation. This ensures that the
        start gain and stop lost mutations before the novel start site are not
        retained. #519 """
        gvf = [
            self.data_dir/'fuzz/09/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/09/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case38(self):
        """ Test case from fuzz test that ensures frameshifting mutations that
        are right after a novel ORF start site and in the same node of it, are
        carried over to downstream nodes. #526 """
        gvf = [
            self.data_dir/'fuzz/10/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/10/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case39(self):
        """ Test case from fuzz test that ensures in-frame deletion, stop
        retaining mutations to be recognized correctly. #527 """
        gvf = [
            self.data_dir/'fuzz/11/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/11/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case40(self):
        """ Test case from fuzz test that ensures peptides with in-frame
        deletion stop lost mutations are called. #528 """
        gvf = [
            self.data_dir/'fuzz/12/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/12/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case41(self):
        """ Test case from fuzz test that ensures variant coordinates being
        handled property when the cleavage site is contained in the inserted
        sequence. #529 """
        gvf = [
            self.data_dir/'fuzz/13/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/13/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case42(self):
        """ Test case from fuzz test with two overlapping deletions. #531 """
        gvf = [
            self.data_dir/'fuzz/14/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/14/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case43(self):
        """ Test case from fuzz test that indel variants after collapsing
        not considered. #533 """
        gvf = [
            self.data_dir/'fuzz/15/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/15/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case44(self):
        """ Test case from fuzz test that nodes are missed when there are
        multiple frameshifting mutation that makes it go back to the original
        reading frame. #534 """
        gvf = [
            self.data_dir/'fuzz/16/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/16/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case45(self):
        """ Test case from fuzz test that stop altering mutation not retained
        during collapsing. #549 """
        gvf = [
            self.data_dir/'fuzz/17/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/17/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case46(self):
        """ Test case from fuzz test that peptides are missed when multiple
        variants cause the same sequence. #552 """
        gvf = [
            self.data_dir/'fuzz/18/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/18/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case47(self):
        """ Test case from fuzz test that `cpop_collapsed` attribute was not
        retained afer merging so peptides that don't end with cleavage sites
        were yield. #552 """
        gvf = [
            self.data_dir/'fuzz/19/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/19/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case48(self):
        """ Test case from fuzz test that subgraph identity got lost after
        mergine. #566 """
        gvf = [
            self.data_dir/'fuzz/20/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/20/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case49(self):
        """ Test case from fuzz test that variants are not filtered correctly
        for fusion. #567 """
        gvf = [
            self.data_dir/'fuzz/21/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/21/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case50(self):
        """ Test case from fuzz test that stop altering variants where there
        is any novel start codon after it should be ignored. #568 """
        gvf = [
            self.data_dir/'fuzz/22/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/22/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case51(self):
        """ An issue caught by fuzz test. In #684b6f we fixe it in the way that
        when finding the farthest downstream reference node for aligning the
        variant bubble, if the downstream node has a incoming node from a different
        subgraph, a farther downstream node will be taken. But in the case of
        fusion with insertions, the first node of a subgraph always have a
        incoming node from a different subgraph. And this causes the entire
        fusion graph to be resolved. So here we made it an exception for fusion.
        """
        gvf = [
            self.data_dir/'fuzz/23/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/23/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case52(self):
        """ Issue found by fuzz test to make sure the correct right most node
        of a variant bubble is found.
        """
        gvf = [
            self.data_dir/'fuzz/24/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/24/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case53(self):
        """ Test case to ensure that because alt splice deletions are not terated
        as subgraphs any more so it won't be returned by `move_downsteams` as
        an end node.
        """
        gvf = [
            self.data_dir/'fuzz/25/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/25/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case54(self):
        """ To ensure that in-frame subgraph won't be treated as subgraph when
        aligning variant bubbles.
        """
        gvf = [
            self.data_dir/'fuzz/26/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/26/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case55(self):
        """ Caught by fuzz test that a start gain mutation is missing in the
        current loop. #576
        """
        gvf = [
            self.data_dir/'fuzz/27/fake_variants.gvf',
            self.data_dir/'fuzz/27/fake_circ_rna.gvf'
        ]
        expected = self.data_dir/'fuzz/27/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_call_variant_peptide_case56(self):
        """ Caught by fuzz test that fusion not inserted into the correct
        position when the breakpoint (eitgher donor or accetper) is intronic.
        #578
        """
        gvf = [
            self.data_dir/'fuzz/28/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/28/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)
