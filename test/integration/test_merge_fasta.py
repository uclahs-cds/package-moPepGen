""" Test mergeFasta """
import argparse
import subprocess as sp
import sys
from test.unit import create_aa_record
from test.integration import TestCaseIntegration
from Bio.SeqIO import FastaIO
from moPepGen import cli
from moPepGen.aa import VariantPeptidePool


class TestMergeFasta(TestCaseIntegration):
    """ Test cases for moPepGen mergeFasta """

    def write_fasta(self, data, file_name):
        """ Write data to fasta file """
        seqs = [create_aa_record(*_) for _ in data]
        record2title = lambda x: x.description
        with open(self.work_dir/file_name, 'wt') as handle:
            writer = FastaIO.FastaWriter(handle, record2title=record2title)
            for seq in seqs:
                writer.write_record(seq)

    def create_base_arg(self):
        """ Create base args """
        args = argparse.Namespace()
        args.command = 'mergeFasta'
        args.output_path = self.work_dir/"output.fasta"
        args.dedup_header = False
        args.denylist = []
        args.quiet = False
        return args

    def test_merge_fasta_cli(self):
        """ Test mergeFasta cli """
        data1 = [
            ('AAAAAK',  'ENST0001|SNV-1-A-T|1'),
            ('AAAAAAK', 'ENST0002|SNV-2-A-T|1')
        ]
        data2 = [
            ('AAAAEK',  'ENST0003|SNV-1-A-T|1'),
            ('AAAAAEK', 'ENST0004|SNV-2-A-T|1')
        ]
        filename1 = 'database1.fasta'
        filename2 = 'database2.fasta'
        self.write_fasta(data1, filename1)
        self.write_fasta(data2, filename2)

        input_path = str(self.work_dir/filename1)
        input_path += f" {self.work_dir/filename2}"

        output_path = str(self.work_dir/'merged.fasta')

        cmd = f"""
        {sys.executable} -m moPepGen.cli mergeFasta -i {input_path} -o {output_path}
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_merge_fasta_case1(self):
        """ Test case when there isn't any conflict """
        data1 = [
            ('AAAAAK',  'ENST0001|SNV-1-A-T|1'),
            ('AAAAAAK', 'ENST0002|SNV-2-A-T|1')
        ]
        data2 = [
            ('AAAAEK',  'ENST0003|SNV-1-A-T|1'),
            ('AAAAAEK', 'ENST0004|SNV-2-A-T|1')
        ]
        filename1 = 'database1.fasta'
        filename2 = 'database2.fasta'
        self.write_fasta(data1, filename1)
        self.write_fasta(data2, filename2)

        args = self.create_base_arg()
        args.input_path = [
            self.work_dir/filename1,
            self.work_dir/filename2
        ]

        cli.merge_fasta(args)

        self.assertTrue(args.output_path)

        with open(args.output_path, 'rt') as handle:
            pool = VariantPeptidePool.load(handle)
        self.assertEqual(len(pool.peptides), 4)

    def test_merge_fasta_case2(self):
        """ Test case when there are conflict of same peptide sequences. """
        data1 = [
            ('AAAAAK',  'ENST0001|SNV-1-A-T|1'),
            ('AAAAAAK', 'ENST0002|SNV-2-A-T|1')
        ]
        data2 = [
            ('AAAAAK',  'ENST0001|SNV-1-A-T|1'),
            ('AAAAAAK', 'ENST0002|SNV-2-A-T|1')
        ]
        filename1 = 'database1.fasta'
        filename2 = 'database2.fasta'
        self.write_fasta(data1, filename1)
        self.write_fasta(data2, filename2)

        args = self.create_base_arg()
        args.input_path = [
            self.work_dir/filename1,
            self.work_dir/filename2
        ]
        args.dedup_header = True

        cli.merge_fasta(args)

        self.assertTrue(args.output_path)

        with open(args.output_path, 'rt') as handle:
            pool = VariantPeptidePool.load(handle)
        self.assertEqual(len(pool.peptides), 2)
        # redundant fasta header entries are removed
        self.assertEqual(
            {x.description for x in pool.peptides},
            {'ENST0001|SNV-1-A-T|1', 'ENST0002|SNV-2-A-T|1'}
        )

    def test_merge_fasta_denylist(self):
        """ Test denylist """
        data1 = [
            ('AAAAAK',  'ENST0001|SNV-1-A-T|1'),
            ('AAAAAAK', 'ENST0002|SNV-2-A-T|1')
        ]
        data2 = [
            ('AAAAAK',  'ENST0001|SNV-1-A-T|1'),
            ('AAAAAAK', 'ENST0002|SNV-2-A-T|1')
        ]
        filename1 = 'database1.fasta'
        filename2 = 'database2.fasta'
        self.write_fasta(data1, filename1)
        self.write_fasta(data2, filename2)

        data3 = [
            ('AAAAAAK', 'ENST0003|ORF1|1')
        ]
        filename3 = 'denylist.fasta'
        self.write_fasta(data3, filename3)

        args = self.create_base_arg()
        args.input_path = [
            self.work_dir/filename1,
            self.work_dir/filename2
        ]
        args.dedup_header = True
        args.denylist = [
            self.work_dir/filename3
        ]

        cli.merge_fasta(args)

        self.assertTrue(args.output_path)

        with open(args.output_path, 'rt') as handle:
            pool = VariantPeptidePool.load(handle)
        self.assertEqual(len(pool.peptides), 1)
        # redundant fasta header entries are removed
        self.assertEqual(
            {str(x.seq) for x in pool.peptides},
            {'AAAAAK'}
        )

    def test_merge_fasta_table(self):
        """ Test merge fasta with peptide table """
        args = self.create_base_arg()
        args.input_path = [
            self.data_dir/'peptides/gsnp.fasta',
            self.data_dir/'peptides/mtsnv.fasta'
        ]

        cli.merge_fasta(args)

        self.assertTrue(args.output_path.exists())
