""" Test the command line interface """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli, seqvar
from moPepGen.cli.common import load_references


class TestParseStarFusion(TestCaseIntegration):
    """ Test cases for moPepGen parseSTARFusion """

    def test_parse_star_fusion_case1(self):
        """ Test parseSTARFusion """
        args = argparse.Namespace()
        args.fusion = self.data_dir/'fusion/star_fusion.txt'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'star_fusion')
        args.verbose = False
        cli.parse_star_fusion(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'star_fusion.tvf'}
        self.assertEqual(files, expected)

        genome, annotation, _ = load_references(args, load_canonical_peptides=False)

        for record in seqvar.io.parse(self.work_dir/'star_fusion.tvf'):
            tx_id = record.location.seqname
            tx_model = annotation.transcripts[tx_id]
            tx_chr = tx_model.transcript.chrom
            tx_seq = tx_model.get_transcript_sequence(genome[tx_chr])
            x = record.location.start
            self.assertEqual(str(tx_seq.seq[x:x+2]), 'GT')

            tx_id = record.attrs['ACCEPTER_TRANSCRIPT_ID']
            tx_model = annotation.transcripts[tx_id]
            tx_chr = tx_model.transcript.chrom
            tx_seq = tx_model.get_transcript_sequence(genome[tx_chr])
            x = record.get_accepter_position()
            self.assertEqual(str(tx_seq.seq[x-2:x]), 'AG')
