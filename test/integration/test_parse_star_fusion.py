""" Test the command line interface """
import argparse
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli, seqvar
from moPepGen.cli.common import load_references


class TestParseStarFusion(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """
    def test_parse_star_fusion_cli(self):
        """ test parseSTARFusion cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli parseSTARFusion \\
            -i {self.data_dir}/fusion/star_fusion.txt \\
            -o {self.work_dir}/fusion_catcher.gvf \\
            -g {self.data_dir}/genome.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            --source Fusion
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_star_fusion_record_case1(self):
        """ Test parseSTARFusion """
        args = argparse.Namespace()
        args.command = 'parseSTARFusion'
        args.source = 'Fusion'
        args.input_path = self.data_dir/'fusion/star_fusion.txt'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.output_path = self.work_dir/'star_fusion.gvf'
        args.min_est_j = 3.0
        args.quiet = True
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'star_fusion.gvf'}
        self.assertEqual(files, expected)

        genome, anno, *_ = load_references(args, load_canonical_peptides=False)

        for record in seqvar.io.parse(self.work_dir/'star_fusion.gvf'):
            gene_id = record.location.seqname
            gene_model = anno.genes[gene_id]
            gene_chr = gene_model.chrom
            gene_seq = gene_model.get_gene_sequence(genome[gene_chr])
            x = record.location.start
            self.assertEqual(str(gene_seq.seq[x:x+2]), 'GT')

            gene_id = record.attrs['ACCEPTER_TRANSCRIPT_ID']
            gene_model = anno.transcripts[gene_id]
            gene_chr = gene_model.transcript.chrom
            gene_seq = gene_model.get_transcript_sequence(genome[gene_chr])
            x = record.get_accepter_position()
            self.assertEqual(str(gene_seq.seq[x-2:x]), 'AG')

    def test_parse_star_fusion_case1(self):
        """ test parseSTARFusion case1 """
        args = argparse.Namespace()
        args.command = 'parseSTARFusion'
        args.input_path = self.data_dir/'fusion/star_fusion.txt'
        args.source = 'Fusion'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.output_path = self.work_dir/'star_fusion.gvf'
        args.min_est_j = 3.0
        args.quiet = True
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'star_fusion.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)
