""" Integration test of splitDatabase """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestSplitDatabase(TestCaseIntegration):
    """ Test cases for splitDatabase """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base args """
        args = argparse.Namespace()
        args.command = 'splitDatabase'
        args.noncoding_peptides = None
        args.order_source = None
        args.group_source = None
        args.max_source_groups = 1
        args.priority_list = None
        args.output_prefix = self.work_dir/'test'
        args.verbose = False
        return args

    def test_split_database_case1(self):
        """ test splitDatabase case1 """
        args = self.create_base_args()
        args.variant_gvf = [
            self.data_dir/'vep/vep_gSNP.gvf',
            self.data_dir/'vep/vep_gINDEL.gvf',
            self.data_dir/'reditools/reditools.gvf',
            self.data_dir/'fusion/star_fusion.gvf'
        ]
        args.variant_peptides = self.data_dir/'peptides/variant.fasta'
        cli.split_database(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'test_gINDEL.fasta','test_gSNP.fasta','test_RNAEditing.fasta'}
        self.assertEqual(files, expected)
