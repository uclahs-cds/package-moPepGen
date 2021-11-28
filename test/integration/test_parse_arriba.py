""" Test parseArriba """
import argparse
from pathlib import Path
from test.unit import load_references
from test.integration import TestCaseIntegration
from moPepGen import cli, parser


class TestParseArriba(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """
    def test_parse_arriba(self):
        """ Test parseFusionCatcher """
        args = argparse.Namespace()
        args.command = 'parseArriba'
        args.fusion = self.data_dir/'fusion/arriba.txt'
        args.source = 'Fusion'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'arriba')
        args.max_common_mapping = 0
        args.min_spanning_unique = 5
        args.verbose = False
        cli.parse_arriba(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'arriba.gvf'}
        self.assertEqual(files, expected)
