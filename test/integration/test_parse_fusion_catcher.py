""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseFusionCatcher(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """
    def test_parse_fusion_catcher(self):
        """ Test parseFusionCatcher """
        args = argparse.Namespace()
        args.fusion = self.data_dir/'fusion/fusion_catcher.txt'
        args.index_dir = self.data_dir/'index'
        args.output_prefix = str(self.work_dir/'fusion_catcher')
        args.verbose = True
        cli.parse_fusion_catcher(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'fusion_catcher.tvf'}
        self.assertEqual(files, expected)
