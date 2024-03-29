""" Test the command line interface """
import argparse
from pathlib import Path
import subprocess as sp
import sys
from test.unit import load_references
from test.integration import TestCaseIntegration
from moPepGen import cli, parser


class TestParseFusionCatcher(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """
    def test_parse_fusioncatcher_cli(self):
        """ test parseFusionCatcher cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli parseFusionCatcher \\
            -i {self.data_dir}/fusion/fusion_catcher.txt \\
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

    def test_fusioncatcher_parser(self):
        """ Test FusionCatcherParser """
        genome, anno = load_references(Path('test/files'))
        parse = parser.FusionCatcherParser.parse
        for record in parse(self.data_dir/'fusion/fusion_catcher.txt'):
            fusion_seq = record.fusion_sequence
            variants = record.convert_to_variant_records(anno, genome)
            for variant in variants:
                gene1_id = variant.location.seqname
                gene1_model = anno.genes[gene1_id]
                chrom1 = gene1_model.chrom
                gene1_seq =  gene1_model.get_gene_sequence(genome[chrom1])
                _end1 = variant.location.start
                _start1 = _end1 - len(fusion_seq[0])
                left_seq = gene1_seq.seq[_start1:_end1]
                self.assertEqual(str(left_seq), fusion_seq[0])

                gene2_id = variant.attrs['ACCEPTER_GENE_ID']
                gene2_model = anno.genes[gene2_id]
                chrom2 = gene2_model.chrom
                gene2_seq =  gene2_model.get_gene_sequence(genome[chrom2])
                _start2 = variant.attrs['ACCEPTER_POSITION']
                _end2 = _start2 + len(fusion_seq[1])
                right_seq = gene2_seq.seq[_start2:_end2]
                self.assertEqual(str(right_seq), fusion_seq[1])

    def test_parse_fusion_catcher(self):
        """ Test parseFusionCatcher """
        args = argparse.Namespace()
        args.command = 'parseFusionCatcher'
        args.input_path = self.data_dir/'fusion/fusion_catcher.txt'
        args.source = 'Fusion'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.reference_source = None
        args.output_path = self.work_dir/'fusion_catcher.gvf'
        args.max_common_mapping = 0
        args.min_spanning_unique = 5
        args.quiet = True
        cli.parse_fusion_catcher(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'fusion_catcher.gvf'}
        self.assertEqual(files, expected)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)
