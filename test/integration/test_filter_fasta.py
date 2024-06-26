""" Integration test module for filterFasta """
import argparse
from pathlib import Path
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from Bio import SeqIO
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
        args.reference_source = None
        args.quiet = True
        args.enzyme = 'trypsin'
        args.miscleavages = None
        args.denylist = None
        args.keep_canonical = False
        return args

    def test_filter_fasta_cli(self):
        """ Test filterFasta cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli filterFasta \\
            -i {self.data_dir}/vep/vep.fasta \\
            -o {self.work_dir}/vep_filtered.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            --exprs-table {self.data_dir}/rsem/rsem.txt \\
            --skip-lines 1 \\
            --tx-id-col 1 \\
            --quant-col 5 \\
            --quant-cutoff 100
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

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

    def test_filter_fasta_miscleavages(self):
        """ test filterFasta to filter miscleavages """
        args = self.generate_default_args()
        args.input_path = Path('test/files/vep/vep.fasta')
        args.exprs_table = None
        args.skip_lines = 0
        args.delimiter = '\t'
        args.tx_id_col = 'transcript_id'
        args.quant_col = 'expected_count'
        args.quant_cutoff = 100
        args.keep_all_coding = False
        args.keep_all_noncoding = False
        args.miscleavages = "2:2"
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)

    def test_filter_fasta_noncoding_case1(self):
        """ test filterFasta to filter novel ORF peptides """
        args = self.generate_default_args()
        args.input_path = Path('test/files/peptides/novel_orf.fasta')
        args.exprs_table = Path('test/files/rsem/rsem.txt')
        args.skip_lines = 1
        args.delimiter = '\t'
        args.tx_id_col = '1'
        args.quant_col = '5'
        args.quant_cutoff = 30
        args.keep_all_coding = False
        args.keep_all_noncoding = False
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)
        peptides_raw = list(SeqIO.parse(args.input_path, 'fasta'))
        peptides_filtered = list(SeqIO.parse(args.output_path, 'fasta'))
        seqs_raw = {x.seq for x in peptides_raw}
        seqs_filtered = {x.seq for x in peptides_filtered}
        self.assertTrue(all(x in seqs_raw for x in seqs_filtered))

    def test_filter_fasta_noncoding_case2(self):
        """ test filterFasta to filer novel ORF peptides that no sequence should
        pass. """
        args = self.generate_default_args()
        args.input_path = Path('test/files/peptides/novel_orf.fasta')
        args.exprs_table = Path('test/files/rsem/rsem.txt')
        args.skip_lines = 1
        args.delimiter = '\t'
        args.tx_id_col = '1'
        args.quant_col = '5'
        args.quant_cutoff = 10000
        args.keep_all_coding = False
        args.keep_all_noncoding = False
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)
        peptides_filtered = list(SeqIO.parse(args.output_path, 'fasta'))
        self.assertEqual(len(peptides_filtered), 0)

    def test_filter_fasta_denylist(self):
        """ Test filterFasta with given denylist """
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
        args.denylist = Path('test/files/peptides/novel_orf.fasta')
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)

    def test_filter_fasta_denylist_only(self):
        """ Test filterFasta with only a given denylist no exprs_table """
        args = self.generate_default_args()
        args.input_path = Path('test/files/vep/vep.fasta')
        args.exprs_table = None
        args.skip_lines = 1
        args.delimiter = '\t'
        args.tx_id_col = '1'
        args.quant_col = '5'
        args.quant_cutoff = 100
        args.keep_all_coding = False
        args.keep_all_noncoding = False
        args.denylist = Path('test/files/peptides/novel_orf.fasta')
        cli.filter_fasta(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'vep_filtered.fasta'}
        self.assertEqual(files, expected)
