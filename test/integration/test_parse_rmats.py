""" Test the moPepGen parseRMATS """
import argparse
from pathlib import Path
import subprocess as sp
import sys
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
        args.reference_source = None
        args.output_path = self.work_dir/'rmats.gvf'
        args.quiet = True
        return args

    def test_parse_rmats_cli(self):
        """ test parseRMATS cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli parseRMATS \\
            --se {self.data_dir}/alternative_splicing/rmats_se_case_1.SE.JC.txt \\
            --a5ss {self.data_dir}/alternative_splicing/rmats_a5ss_case_1.A5SS.JC.txt \\
            --a3ss {self.data_dir}/alternative_splicing/rmats_a3ss_case_1.A3SS.JC.txt \\
            --mxe {self.data_dir}/alternative_splicing/rmats_mxe_case_1.MXE.JC.txt \\
            --ri {self.data_dir}/alternative_splicing/rmats_ri_case_1.RI.JC.txt \\
            -o {self.work_dir}/reditools.gvf \\
            -g {self.data_dir}/genome.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            --source rMATS
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_parse_rmats_se_case_1(self):
        """ rMATS skipped exon when the retained version is annotated. This
        should results an deletion. """
        args = self.create_base_args()
        args.skipped_exon = self.data_dir\
            /'alternative_splicing/rmats_se_case_1.SE.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.location.start, 323)
        self.assertTrue(record.id, 'SE_324')
        self.assertTrue(int(record.attrs['START']), 323)
        self.assertTrue(int(record.attrs['END']), 405)
        self.assertTrue(record.type, 'Deletion')
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_rmats_se_case_2(self):
        """ rMATS skipped exon when the skipped version is annotated. This
        should result an insertion. """
        args = self.create_base_args()
        args.skipped_exon = self.data_dir\
            /'alternative_splicing/rmats_se_case_2.SE.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.location.start, 870)
        self.assertTrue(record.id, 'SE_870')
        self.assertTrue(int(record.attrs['DONOR_START']), 870)
        self.assertTrue(int(record.attrs['DONOR_END']), 1097)
        self.assertTrue(record.type, 'Insertion')
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_rmats_se_case_3(self):
        """ Complex cases with:
        1. Multiple exons skipped.
        2. Upstream adjacent exons
        3. Downstream adjacent exons
        """
        args = self.create_base_args()
        args.skipped_exon = self.data_dir\
            /'alternative_splicing/rmats_se_case_3.SE.JC.txt'
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(args.output_path))
        expect_values = {
            ('Deletion', 323, 405),
            ('Deletion', 323, 750),
            ('Deletion', 405, 869),
            ('Deletion', 750, 869),
            ('Substitution', 750, 1097, 1097, 1264),
            ('Substitution', 323, 405, 405, 750),
            ('Insertion', 404, 405, 750)
        }
        received_values = set()
        for record in records:
            if record.type == 'Deletion':
                received_values.add(
                    ('Deletion', int(record.attrs['START']), int(record.attrs['END']))
                )
            elif record.type == 'Substitution':
                received_values.add((
                    'Substitution', int(record.attrs['START']), int(record.attrs['END']),
                    int(record.attrs['DONOR_START']), int(record.attrs['DONOR_END'])
                ))
            else:
                received_values.add((
                    'Insertion', record.location.start,
                    int(record.attrs['DONOR_START']), int(record.attrs['DONOR_END'])
                ))
        self.assertEqual(received_values, expect_values)

    def test_parse_rmats_a5ss_case_1(self):
        """ rMATS A5SS when the longer version is annotated. This should
        results a deletion. """
        args = self.create_base_args()
        args.alternative_5_splicing = self.data_dir\
            /'alternative_splicing/rmats_a5ss_case_1.A5SS.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.type, 'Deletion')
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_rmats_a5ss_case_2(self):
        """ rMATS A5SS when the shorter version is annotated. This should
        results an insertion. """
        args = self.create_base_args()
        args.alternative_5_splicing = self.data_dir\
            /'alternative_splicing/rmats_a5ss_case_2.A5SS.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.type, 'Insertion')
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_rmats_a3ss_case_1(self):
        """ rMATS A3SS when the longer version is annotated. This should
        results a deletion. """
        args = self.create_base_args()
        args.alternative_3_splicing = self.data_dir\
            /'alternative_splicing/rmats_a3ss_case_1.A3SS.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.type, 'Deletion')
        self.assert_gvf_order(args.output_path, args.annotation_gtf)

    def test_parse_rmats_a3ss_case_2(self):
        """ rMATS A3SS when the shorter version is annotated. This should
        results an Insertion. """
        args = self.create_base_args()
        args.alternative_3_splicing = self.data_dir\
            /'alternative_splicing/rmats_a3ss_case_2.A3SS.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.type, 'Insertion')
        self.assert_gvf_order(Path(args.output_path), args.annotation_gtf)

    def test_parse_rmats_mxe_case_1(self):
        """ rMATS MXE when one exon is annotated. This should results a
        substitution. """
        args = self.create_base_args()
        args.mutually_exclusive_exons = self.data_dir\
            /'alternative_splicing/rmats_mxe_case_1.MXE.JC.txt'
        cli.parse_rmats(args)
        record = list(seqvar.io.parse(args.output_path))[0]
        self.assertTrue(record.type, 'Substitution')
        self.assert_gvf_order(Path(args.output_path), args.annotation_gtf)

    def test_parse_rmats_mxe_case_2(self):
        """ rMATS MXE when both exons are annotated. Should be two deletions """
        args = self.create_base_args()
        args.mutually_exclusive_exons = self.data_dir\
            /'alternative_splicing/rmats_mxe_case_2.MXE.JC.txt'
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(args.output_path))
        self.assertTrue(all(v.type == 'Deletion' for v in records))

    def test_parse_rmats_ri(self):
        """ rMATS RI. """
        args = self.create_base_args()
        args.retained_intron = self.data_dir\
            /'alternative_splicing/rmats_ri_case_1.RI.JC.txt'
        cli.parse_rmats(args)
        records = list(seqvar.io.parse(args.output_path))
        self.assertTrue(len(records), 2)
        self.assert_gvf_order(args.output_path, args.annotation_gtf)
        for record in records:
            self.assertEqual(record.type, 'Insertion')
