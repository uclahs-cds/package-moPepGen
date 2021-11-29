""" Test parseArriba """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


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
        args.min_split_read1 = 1
        args.min_split_read2 = 1
        args.min_confidence = 'medium'
        args.verbose = False
        cli.parse_arriba(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'arriba.gvf'}
        self.assertEqual(files, expected)
