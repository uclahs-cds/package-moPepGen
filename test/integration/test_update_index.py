""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestGenerateIndex(TestCaseIntegration):
    """ Test cases for moPepGen generateIndex """
    def setUp(self):
        super().setUp()
        args = argparse.Namespace()
        args.command = 'generateIndex'
        args.genome_fasta = self.data_dir / 'genome.fasta'
        args.annotation_gtf = self.data_dir / 'annotation.gtf'
        args.proteome_fasta = self.data_dir / 'translate.fasta'
        args.codon_table = 'Standard'
        args.chr_codon_table = ['chrM:SGC1']
        args.start_codons = ['ATG']
        args.chr_start_codons = ['chrM:ATG,ATA,ATT']
        args.gtf_symlink = False
        args.reference_source = None
        args.invalid_protein_as_noncoding = False
        args.cleavage_rule = 'trypsin'
        args.cleavage_exception = 'trypsin_exception'
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.miscleavage = 2
        args.quiet = True
        args.force = False
        args.output_dir = self.work_dir / 'index'
        cli.generate_index(args)

    def create_base_args(self) -> argparse.Namespace:
        """ create base args """
        args = argparse.Namespace()
        args.command = 'updateIndex'
        args.gtf_symlink = False
        args.reference_source = None
        args.invalid_protein_as_noncoding = False
        args.cleavage_rule = 'trypsin'
        args.cleavage_exception = 'trypsin_exception'
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.miscleavage = 2
        args.quiet = True
        args.force = False
        args.index_dir = self.work_dir / 'index'
        return args

    def test_update_index_case1(self):
        """ test udpate index case 1 """
        args = self.create_base_args()
        args.cleavage_rule = 'lysc'
        cli.update_index(args)

    def test_update_index_case2(self):
        """ test udpate index case 2 """
        args = self.create_base_args()
        args.cleavage_rule = 'trypsin'
        with self.assertRaises(SystemExit):
            cli.update_index(args)

    def test_update_index_case3(self):
        """ test udpate index case 1 """
        args = self.create_base_args()
        args.cleavage_rule = 'trypsin'
        args.force = True
        cli.update_index(args)
        self.assertEqual(
            {x.name for x in args.index_dir.rglob('canonical_peptides*')},
            {'canonical_peptides_001.pkl'}
        )

    def test_update_index_no_cleavage(self):
        """ test udpate index case 1 """
        args = self.create_base_args()
        args.cleavage_rule = None
        cli.update_index(args)
        # No additional canonical peptides file should be created
        self.assertEqual(
            {x.name for x in args.index_dir.rglob('canonical_peptides*')},
            {'canonical_peptides_001.pkl'}
        )
