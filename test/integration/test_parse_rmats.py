""" Test the moPepGen parseRMATS """
import argparse
from test.integration import TestCaseIntegration
from moPepGen import cli, seqvar


class TestParseRMATS(TestCaseIntegration):
    """ Test cases for moPepGen parseRMATS """

    def create_base_args(self) -> argparse.Namespace:
        """ Create base args """
        args = argparse.Namespace()
        args.command = 'parseRMATS'
        args.source = 'AlternativeSplicing'
        args.skipped_exon = None
        args.alternative_5_splicing = None
        args.alternative_3_splicing = None
        args.mutually_exclusive_exons = None
        args.retained_intron = None
        args.min_sjc = 1
        args.min_ijc = 1
        args.index_dir = None
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_prefix = str(self.work_dir/'rmats')
        args.verbose = False
        return args

    def test_parse_rmats_se_case_1(self):
        """ rMATS skipped exon when the retained version is annotated. This
        should results an deletion. """
        args = self.create_base_args()
        args.skipped_exon = self.data_dir/'alternative_splicing/rmats_se_case_1.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.location.start, 323)
        self.assertTrue(record.id, 'SE_324')
        self.assertTrue(int(record.attrs['START']), 323)
        self.assertTrue(int(record.attrs['END']), 405)
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_se_case_2(self):
        """ rMATS skipped exon when the skipped version is annotated. This
        should result an insertion. """
        args = self.create_base_args()
        args.skipped_exon = self.data_dir/'alternative_splicing/rmats_se_case_2.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.location.start, 870)
        self.assertTrue(record.id, 'SE_870')
        self.assertTrue(int(record.attrs['DONOR_START']), 870)
        self.assertTrue(int(record.attrs['DONOR_END']), 1097)
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_a5ss_case_1(self):
        """ rMATS A5SS when the longer version is annotated. This should
        results a deletion. """
        args = self.create_base_args()
        args.alternative_5_splicing = self.data_dir/'alternative_splicing/rmats_a5ss_case_1.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_a5ss_case_2(self):
        """ rMATS A5SS when the shorter version is annotated. This should
        results an insertion. """
        args = self.create_base_args()
        args.alternative_5_splicing = self.data_dir/'alternative_splicing/rmats_a5ss_case_2.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_a3ss_case_1(self):
        """ rMATS A3SS when the longer version is annotated. This should
        results a deletion. """
        args = self.create_base_args()
        args.alternative_3_splicing = self.data_dir/'alternative_splicing/rmats_a3ss_case_1.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.type, 'Deletion')

    def test_parse_rmats_a3ss_case_2(self):
        """ rMATS A3SS when the shorter version is annotated. This should
        results an Insertion. """
        args = self.create_base_args()
        args.alternative_3_splicing = self.data_dir/'alternative_splicing/rmats_a3ss_case_2.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.type, 'Insertion')

    def test_parse_rmats_mxe_case_1(self):
        """ rMATS MXE when one exon is annotated. This should results a
        substitution. """
        args = self.create_base_args()
        args.mutually_exclusive_exons = self.data_dir/'alternative_splicing/rmats_mxe_case_1.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))[0]
        self.assertTrue(record.type, 'Substitution')

    def test_parse_rmats_mxe_case_2(self):
        """ rMATS MXE when both exons are annotated. This should results two
        deletions. """
        args = self.create_base_args()
        args.mutually_exclusive_exons = self.data_dir/'alternative_splicing/rmats_mxe_case_2.txt'
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))
        self.assertTrue(len(records), 2)
        for record in records:
            self.assertEqual(record.type, 'Deletion')

    def test_parse_rmats_ri(self):
        """ rMATS RI. """
        args = self.create_base_args()
        args.retained_intron = self.data_dir/'alternative_splicing/rmats_ri_case_1.txt'
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(f'{args.output_prefix}.gvf'))
        self.assertTrue(len(records), 2)
        for record in records:
            self.assertEqual(record.type, 'Insertion')
