""" Test module for parseREDItools """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseREDItools(TestCaseIntegration):
    """ Test parseREDItools """
    def test_parse_reditools_case1(self):
        """ Test parse reditools """
        args = argparse.Namespace()
        args.command = 'parseREDItools'
        args.source = 'RNAEditingSite'
        args.reditools_table = self.data_dir/'reditools/reditools_annotated.txt'
        args.transcript_id_column = 16
        args.index_dir = None
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'reditools')
        args.quiet = True
        args.min_read_count = 3
        args.min_frequency = 0.1
        cli.parse_reditools(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'reditools.gvf'}
        self.assertEqual(files, expected)
