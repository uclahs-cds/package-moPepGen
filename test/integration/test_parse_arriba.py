""" Test parseArriba """
import argparse
import subprocess as sp
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseArriba(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """
    def test_parse_arriba(self):
        """ Test parseArriba """
        args = argparse.Namespace()
        args.command = 'parseArriba'
        args.input_path = self.data_dir/'fusion/arriba.txt'
        args.source = 'Fusion'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'arriba.gvf'
        args.min_split_read1 = 1
        args.min_split_read2 = 1
        args.min_confidence = 'medium'
        args.quiet = False
        cli.parse_arriba(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'arriba.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_arriba_cli(self):
        """ Test parseArriba from command line """
        cmd = f"""
        python -m moPepGen.cli parseArriba \
            -i {self.data_dir}/fusion/arriba.txt \
            -o {self.work_dir}/arriba.gvf \
            -g {self.data_dir}/genome.fasta \
            -a {self.data_dir}/annotation.gtf \
            -p {self.data_dir}/translate.fasta
        """
        sp.run(cmd, shell=True, check=False)
