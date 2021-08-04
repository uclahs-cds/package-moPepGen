""" Test the moPepGen parseRMATS """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli, seqvar


class TestParseRMATS(TestCaseIntegration):
    """ Test cases for moPepGen parseRMATS """

    def test_parse_rmats_se_case_1(self):
        """ rMATS skipped exon when the retained version is annotated. This
        should results an deletion. """
        args = argparse.Namespace()
        args.skipped_exon = self.data_dir/'alternative_splicing/rmats_se_case_1.txt'
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.location.start, 323)
        self.assertTrue(record.id, 'SE_324')
        self.assertTrue(int(record.attrs['START']), 323)
        self.assertTrue(int(record.attrs['END']), 405)
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_se_case_2(self):
        """ rMATS skipped exon when the skipped version is annotated. This
        should result an insertion. """
        args = argparse.Namespace()
        args.skipped_exon = self.data_dir/'alternative_splicing/rmats_se_case_2.txt'
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.location.start, 870)
        self.assertTrue(record.id, 'SE_870')
        self.assertTrue(int(record.attrs['START']), 870)
        self.assertTrue(int(record.attrs['END']), 1097)
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_a5ss_case_1(self):
        """ rMATS A5SS when the longer version is annotated. This should
        results a deletion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = self.data_dir/'alternative_splicing/rmats_a5ss_case_1.txt'
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_a5ss_case_2(self):
        """ rMATS A5SS when the shorter version is annotated. This should
        results an insertion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = self.data_dir/'alternative_splicing/rmats_a5ss_case_2.txt'
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_a3ss_case_1(self):
        """ rMATS A3SS when the longer version is annotated. This should
        results a deletion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = self.data_dir/'alternative_splicing/rmats_a3ss_case_1.txt'
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_a3ss_case_2(self):
        """ rMATS A3SS when the shorter version is annotated. This should
        results an Insertion. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = self.data_dir/'alternative_splicing/rmats_a3ss_case_2.txt'
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_mxe_case_1(self):
        """ rMATS MXE when one exon is annotated. This should results a
        substitution. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = self.data_dir/'alternative_splicing/rmats_mxe_case_1.txt'
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))[0]
        self.assertTrue(record.type, 'Substitution')

    def test_parse_rmats_mxe_case_2(self):
        """ rMATS MXE when both exons are annotated. This should results two
        deletions. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = self.data_dir/'alternative_splicing/rmats_mxe_case_2.txt'
        args.retained_intron = None
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))
        self.assertTrue(len(records), 2)
        for record in records:
            self.assertEqual(record.type, 'Deletion')

    def test_parse_rmats_ri(self):
        """ rMATS RI. """
        args = argparse.Namespace()
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = self.data_dir/'alternative_splicing/rmats_ri_case_1.txt'
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(f'{args.output_prefix}.tvf'))
        self.assertTrue(len(records), 2)
        for record in records:
            self.assertEqual(record.type, 'Insertion')
