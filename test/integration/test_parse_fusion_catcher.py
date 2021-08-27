""" Test the command line interface """
import argparse
from pathlib import Path
from test.unit import load_references
from test.integration import TestCaseIntegration
from moPepGen import cli, parser


class TestParseFusionCatcher(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """
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
        args.fusion = self.data_dir/'fusion/fusion_catcher.txt'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.output_prefix = str(self.work_dir/'fusion_catcher')
        args.max_common_mapping = 0
        args.min_spanning_unique = 5
        args.verbose = False
        cli.parse_fusion_catcher(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'fusion_catcher.gvf'}
        self.assertEqual(files, expected)
