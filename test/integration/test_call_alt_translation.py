""" Test the CLI for callAltTranslation """
import argparse
import subprocess as sp
import sys
import functools
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli
from moPepGen.aa import VariantPeptideIdentifier as vpi


def create_base_args() -> argparse.Namespace:
    """ Create a base args """
    args = argparse.Namespace()
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.reference_source = None
    args.output_path = None
    args.w2f_reassignment = True
    args.selenocysteine_termination = True
    args.cleavage_rule = 'trypsin'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.quiet = True
    return args

class TestCallAltTranslation(TestCaseIntegration):
    """ Test cases for moPepGen callAltTranslation """

    def test_call_alt_translation_cli(self):
        """ test callAltTranslation cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli callAltTranslation \\
            -o {self.work_dir}/alt_translation.fasta \\
            -g {self.data_dir}/genome.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta \\
            --w2f-reassignment \\
            --selenocysteine-termination
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_call_alt_translation_case1(self):
        """ test call alt translation peptides """
        args = create_base_args()
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'alt_translation.fasta'

        cli.call_alt_translation(args)
        with open(args.output_path, 'rt') as handle:
            peptides = list(SeqIO.parse(handle, 'fasta'))
            headers = [p.description for p in peptides]
            self.assertEqual(len(headers),len(set(headers)))

            var_labels = functools.reduce(
                lambda x,y: x + y,
                [vpi.parse_variant_peptide_id(x) for x in headers]
            )
            self.assertTrue(all(isinstance(x, vpi.BaseVariantPeptideIdentifier)
                for x in var_labels))
            self.assertTrue(all(len(x.variant_ids) > 0 for x in var_labels))
