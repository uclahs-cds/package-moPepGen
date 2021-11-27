""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestGenerateIndex(TestCaseIntegration):
    """ Test cases for moPepGen generateIndex """

    def test_generate_index_case1(self):
        """ Test genreate index """
        args = argparse.Namespace()
        args.genome_fasta = self.data_dir / 'genome.fasta'
        args.annotation_gtf = self.data_dir / 'annotation.gtf'
        args.proteome_fasta = self.data_dir / 'translate.fasta'
        args.cleavage_rule = 'trypsin'
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.miscleavage = 2
        args.quiet = True
        args.output_dir = self.work_dir / 'index'
        args.output_dir.mkdir(parents=False, exist_ok=True)
        cli.generate_index(args)
        files = {str(file.name) for file in args.output_dir.glob('*.pickle')}
        expected = {'genome.pickle', 'proteome.pickle', 'annotation.pickle',
            'canonical_peptides.pickle'}
        self.assertEqual(files, expected)
