""" Test moPepGen-util brueForce"""
import argparse
from typing import List
import io
import contextlib
from pathlib import Path
from test.integration import TestCaseIntegration
from moPepGen import util


def create_base_args() -> argparse.Namespace:
    """ Create base args """
    args = argparse.Namespace()
    args.cleavage_rule = 'trypsin'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.variant_ids = []
    args.force = True
    return args

class TestBruteForce(TestCaseIntegration):
    """ Test cases for moPepGen-util bruteForce """

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
        args.input_gvf = gvf
        args.reference_dir = reference
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            util.brute_force(args)
            seqs = set(stream.getvalue().rstrip().split('\n'))

        with open(expect, 'rt') as handle:
            expected = {line.strip() for line in handle}

        self.assertEqual(seqs, expected)

    def test_brute_force_fuzz_1(self):
        """ Fuzz test 1 """
        gvf = [
            self.data_dir/'fuzz/01/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/01/brute_force.fasta'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_2(self):
        """ Fuzz test 2 """
        gvf = [
            self.data_dir/'fuzz/02/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/02/brute_force.fasta'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_3(self):
        """ Fuzz test 3 """
        gvf = [
            self.data_dir/'fuzz/03/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/03/brute_force.fasta'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_4(self):
        """ Fuzz test 4 """
        gvf = [
            self.data_dir/'fuzz/04/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/04/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_5(self):
        """ Fuzz test 5 """
        gvf = [
            self.data_dir/'fuzz/05/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/05/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_6(self):
        """ Fuzz test 6 """
        gvf = [
            self.data_dir/'fuzz/06/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/06/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_7(self):
        """ Fuzz test 7 """
        gvf = [
            self.data_dir/'fuzz/07/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/07/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_8(self):
        """ Fuzz test 8 """
        gvf = [
            self.data_dir/'fuzz/08/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/08/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_9(self):
        """ Fuzz test 9 """
        gvf = [
            self.data_dir/'fuzz/09/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/09/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_10(self):
        """ Fuzz test 10 """
        gvf = [
            self.data_dir/'fuzz/10/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/10/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_11(self):
        """ Fuzz test 11 """
        gvf = [
            self.data_dir/'fuzz/11/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/11/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_12(self):
        """ Fuzz test 12 """
        gvf = [
            self.data_dir/'fuzz/12/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/12/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_13(self):
        """ Fuzz test 13 """
        gvf = [
            self.data_dir/'fuzz/13/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/13/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_14(self):
        """ Fuzz test 14 """
        gvf = [
            self.data_dir/'fuzz/14/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/14/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_15(self):
        """ Fuzz test 15 """
        gvf = [
            self.data_dir/'fuzz/15/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/15/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_16(self):
        """ Fuzz test 16 """
        gvf = [
            self.data_dir/'fuzz/16/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/16/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000314675.11'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_17(self):
        """ Fuzz test 17 """
        gvf = [
            self.data_dir/'fuzz/17/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/17/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_18(self):
        """ Fuzz test 18 """
        gvf = [
            self.data_dir/'fuzz/18/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/18/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_19(self):
        """ Fuzz test 19 """
        gvf = [
            self.data_dir/'fuzz/19/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/19/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000452737.5'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_20(self):
        """ Fuzz test 20 """
        gvf = [
            self.data_dir/'fuzz/20/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/20/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_21(self):
        """ Fuzz test 21 """
        gvf = [
            self.data_dir/'fuzz/21/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/21/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_22(self):
        """ Fuzz test 22 """
        gvf = [
            self.data_dir/'fuzz/22/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/22/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_23(self):
        """ Fuzz test 23 """
        gvf = [
            self.data_dir/'fuzz/23/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/23/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_24(self):
        """ Fuzz test 24 """
        gvf = [
            self.data_dir/'fuzz/24/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/24/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_25(self):
        """ Fuzz test 25 """
        gvf = [
            self.data_dir/'fuzz/25/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/25/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_26(self):
        """ Fuzz test 26 """
        gvf = [
            self.data_dir/'fuzz/26/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/26/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_27(self):
        """ Fuzz test 27 """
        gvf = [
            self.data_dir/'fuzz/27/fake_variants.gvf',
            self.data_dir/'fuzz/27/fake_circ_rna.gvf'
        ]
        expected = self.data_dir/'fuzz/27/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_28(self):
        """ Fuzz test 28 """
        gvf = [
            self.data_dir/'fuzz/28/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/28/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_29(self):
        """ Fuzz test 29 """
        gvf = [
            self.data_dir/'fuzz/29/fake_variants.gvf'
        ]
        expected = self.data_dir/'fuzz/29/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)

    def test_brute_force_fuzz_30(self):
        """ Fuzz test 30 """
        gvf = [
            self.data_dir/'fuzz/30/fake_variants.gvf',
            self.data_dir/'fuzz/30/fake_circ_rna.gvf'
        ]
        expected = self.data_dir/'fuzz/30/brute_force.txt'
        reference = self.data_dir/'downsampled_reference/ENST00000265138.4-ENST00000650150.1'
        self.default_test_case(gvf, reference, expected)
