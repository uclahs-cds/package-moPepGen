""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseStarFusion(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """

    def test_parse_star_fusion_case1(self):
        """ Test parseSTARFusion """
        args = argparse.Namespace()
        args.fusion = self.data_dir/'fusion/star_fusion.txt'
        args.index_dir = self.data_dir/'index'
        args.output_prefix = str(self.work_dir/'star_fusion')
        args.verbose = True
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'star_fusion.tvf'}
        self.assertEqual(files, expected)
