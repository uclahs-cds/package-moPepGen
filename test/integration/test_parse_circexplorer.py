""" Integration test for parseCIRCexplorer """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseCIRCexplorer(TestCaseIntegration):
    """ Test cases for moPepGen parseCIRCexplorer """
    def test_parse_circexplorer2(self):
        """ Test parseCIRCexplorer """
        args = argparse.Namespace()
        args.command = 'parseCIRCexplorer'
        args.input_path = self.data_dir/'circRNA/CIRCexplorer_circularRNA_known.txt'
        args.output_prefix = str(self.work_dir/'circ')
        args.source = 'circRNA'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.circexplorer3 = False
        args.min_read_number = 1
        args.verbose = False
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.gvf'}
        self.assertEqual(files, expected)

    def test_parse_circexplorer3(self):
        """ Test parseCIRCexplorer for CIRCexplorer3 """
        args = argparse.Namespace()
        args.command = 'parseCIRCexplorer'
        args.input_path = self.data_dir/'circRNA/CIRCexplorer3_circularRNA_known.txt'
        args.output_prefix = str(self.work_dir/'circ')
        args.source = 'circRNA'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.circexplorer3 = True
        args.min_read_number = 1
        args.min_fbr_circ = None
        args.min_circ_score = None
        args.verbose = False
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.gvf'}
        self.assertEqual(files, expected)

        args.min_fbr_circ = 1
        args.min_circ_score = 1
        cli.parse_circexplorer(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'circ.gvf'}
        self.assertEqual(files, expected)
