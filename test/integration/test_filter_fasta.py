""" Integration test module for filterFasta """
import argparse
from pathlib import Path
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestFilterFasta(TestCaseIntegration):
    """ Test cases for filterFasta """
    def generate_default_args(self) -> argparse.Namespace:
        """ Generate default args """
        args = argparse.Namespace()
        args.command = 'fitlerFasta'
        args.output_path = self.work_dir/'vep_filtered.fasta'
        args.index_dir = None
        args.genome_fasta = Path('test/files/genome.fasta')
        args.annotation_gtf = Path('test/files/annotation.gtf')
        args.proteome_fasta = Path('test/files/translate.fasta')
        args.quiet = True
        return args

    def test_filter_fasta_int_index(self):
        """ test filterFasta """
        args = self.generate_default_args()
        args.input_path = Path('test/files/vep/vep.fasta')
        args.exprs_table = Path('test/files/rsem/rsem.txt')
        args.skip_lines = 1
        args.delimiter = '\t'
        args.tx_id_col = '1'
        args.quant_col = '5'
        args.quant_cutoff = 100
        args.keep_all_coding = False
        args.keep_all_noncoding = False
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)

    def test_filter_fasta_str_index(self):
        """ test filterFasta """
        args = self.generate_default_args()
        args.input_path = Path('test/files/vep/vep.fasta')
        args.exprs_table = Path('test/files/rsem/rsem.txt')
        args.skip_lines = 0
        args.delimiter = '\t'
        args.tx_id_col = 'transcript_id'
        args.quant_col = 'expected_count'
        args.quant_cutoff = 100
        args.keep_all_coding = False
        args.keep_all_noncoding = False
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)
