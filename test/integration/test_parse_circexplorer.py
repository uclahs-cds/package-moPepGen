""" Integration test for parseCIRCexplorer """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseCIRCexplorer(TestCaseIntegration):
    """ Test cases for moPepGen parseCIRCexplorer """
    def test_parse_circexplorer(self):
        """ Test parseCIRCexplorer """
        args = argparse.Namespace()
        args.input_path = self.data_dir/'circRNA/CIRCexplorer_circularRNA_known.txt'
        args.output_prefix = str(self.work_dir/'circ')
        args.source = 'circRNA'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.min_read_number = 1
        args.verbose = False
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.tsv'}
        self.assertEqual(files, expected)
