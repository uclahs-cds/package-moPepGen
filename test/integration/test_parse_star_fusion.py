""" Test the command line interface """
import argparse
import subprocess as sp
import sys
from unittest import mock
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
            -p {self.data_dir}/translate.fasta \\
            --source Fusion
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def create_base_args(self) -> argparse.Namespace:
        """ Create base args """
        args = argparse.Namespace()
        args.command = 'parseSTARFusion'
        args.source = 'Fusion'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.reference_source = None
        args.output_path = self.work_dir/'star_fusion.gvf'
        args.min_est_j = 3.0
        args.quiet = True
        args.skip_failed =False
        return args

    def test_star_fusion_record_case1(self):
        """ Test parseSTARFusion """
        args = self.create_base_args()
        args.input_path = self.data_dir/'fusion/star_fusion.txt'
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'star_fusion.gvf'}
        self.assertEqual(files, expected)

        ref_data = load_references(args, load_canonical_peptides=False)
        genome = ref_data.genome
        anno = ref_data.anno

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
        args = self.create_base_args()
        args.input_path = self.data_dir/'fusion/star_fusion.txt'
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'star_fusion.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    @mock.patch(
        "moPepGen.parser.STARFusionParser.STARFusionRecord.convert_to_variant_records",
        new=mock.MagicMock(side_effect=ValueError())
    )
    def test_parse_star_fusion_skip_failed(self):
        """ test parseSTARFusion case1 """
        args = self.create_base_args()
        args.input_path = self.data_dir/'fusion/star_fusion.txt'
        with self.assertRaises(ValueError):
            cli.parse_star_fusion(args)
        args.skip_failed = True
        cli.parse_star_fusion(args)
