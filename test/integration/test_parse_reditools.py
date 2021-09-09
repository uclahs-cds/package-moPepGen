""" Test module for parseREDItools """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseREDItools(TestCaseIntegration):
    """ Test parseREDItools """
    def test_parse_reditools_case1(self):
        """ Test parse reditools """
        args = argparse.Namespace()
        args.source = 'RNAEditingSite'
        args.reditools_table = self.data_dir/'reditools/reditools_annotated.txt'
        args.transcript_id_column = 16
        args.index_dir = None
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.output_prefix = str(self.work_dir/'reditools')
        args.verbose = False
        cli.parse_reditools(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'reditools.gvf'}
        self.assertEqual(files, expected)
