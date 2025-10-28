""" Test the CLI for callNovelORF """
import argparse
import os
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from Bio import SeqIO
from moPepGen import cli, gtf, dna, aa


def create_base_args() -> argparse.Namespace:
    """ Create a base args """
    args = argparse.Namespace()
    args.command = 'callNovelORF'
    args.index_dir = None
    args.genome_fasta = None
    args.annotation_gtf = None
    args.proteome_fasta = None
    args.reference_source = None
    args.codon_table = 'Standard'
    args.chr_codon_table = ['chrM:SGC1']
    args.start_codons = ['ATG']
    args.chr_start_codons = ['chrM:ATG,ATA,ATT']
    args.output_path = None
    args.output_orf = None
    args.coding_novel_orf = False
    args.inclusion_biotypes = None
    args.exclusion_biotypes = None
    args.coding_novel_orf = False
    args.min_tx_length = 21
    args.orf_assignment = 'max'
    args.w2f_reassignment = False
    args.cleavage_rule = 'trypsin'
    args.cleavage_exception = 'trypsin_exception'
    args.miscleavage = '2'
    args.min_mw = '500.'
    args.min_length = 7
    args.max_length = 25
    args.quiet = True
    return args

class TestCallNovelORFPeptides(TestCaseIntegration):
    """ Test cases for moPepGen callNovelORF """

    def test_call_novel_orf_cli(self):
        """ test callNovelORF cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli callNovelORF \\
            -o {self.work_dir}/circ.fasta \\
            -g {self.data_dir}/genome.fasta \\
            -a {self.data_dir}/annotation.gtf \\
            -p {self.data_dir}/translate.fasta
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_call_novel_orf_peptides_case1(self):
        """ test call novel ORF peptides """
        args = create_base_args()
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        cli.call_novel_orf_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding_peptide.fasta', 'noncoding_orf.fasta'}
        self.assertEqual(files, expected)
        with open(self.work_dir/'noncoding_peptide.fasta', 'rt') as handle:
            peptides = list(SeqIO.parse(handle, 'fasta'))
            ids = [p.id for p in peptides]
            self.assertEqual(len(ids),len(set(ids)))
            self.assertTrue(peptides[0].id.split('|')[0] != 'None')

    def test_call_novel_orf_peptides_case2(self):
        """ test call novel ORF peptides when no ORF is found """
        args = create_base_args()
        ref_dir = self.data_dir/'downsampled_reference/ENST00000644482.1'
        args.genome_fasta = ref_dir/'genome.fasta'
        args.annotation_gtf = ref_dir/'annotation.gtf'
        args.proteome_fasta = ref_dir/'proteome.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        cli.call_novel_orf_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding_peptide.fasta', 'noncoding_orf.fasta'}
        self.assertEqual(files, expected)
        size = os.stat(self.work_dir/'noncoding_peptide.fasta').st_size
        self.assertEqual(size, 0)
        size = os.stat(self.work_dir/'noncoding_orf.fasta').st_size
        self.assertEqual(size, 0)

    def test_call_novel_orf_peptides_w2f(self):
        """ With w2f reassignment """
        args = create_base_args()
        args.w2f_reassignment = True
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        cli.call_novel_orf_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {'noncoding_peptide.fasta', 'noncoding_orf.fasta'}
        self.assertEqual(files, expected)
        with open(self.work_dir/'noncoding_peptide.fasta', 'rt') as handle:
            peptides = list(SeqIO.parse(handle, 'fasta'))
            ids = [p.id for p in peptides]
            self.assertEqual(len(ids),len(set(ids)))
            self.assertTrue(peptides[0].id.split('|')[0] != 'None')
            self.assertTrue(any('W2F' in x for x in ids))

    def test_call_novel_orf_peptides_coding(self):
        """ Test calling for novel ORF peptides from coding TXs """
        args = create_base_args()
        args.w2f_reassignment = True
        args.genome_fasta = self.data_dir/'genome.fasta'
        args.annotation_gtf = self.data_dir/'annotation.gtf'
        args.proteome_fasta = self.data_dir/'translate.fasta'
        args.output_path = self.work_dir/'noncoding_peptide.fasta'
        args.output_orf = self.work_dir/'noncoding_orf.fasta'
        args.coding_novel_orf = True
        cli.call_novel_orf_peptide(args)
        files = {str(file.name) for file in self.work_dir.glob('*')}
        expected = {args.output_path.name, args.output_orf.name}
        self.assertEqual(files, expected)

        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(args.annotation_gtf)

        genome = dna.DNASeqDict()
        genome.dump_fasta(args.genome_fasta)

        proteome = aa.AminoAcidSeqDict()
        proteome.dump_fasta(args.proteome_fasta)
        canonical_peptides = proteome.create_unique_peptide_pool(
            anno=anno, rule='trypsin', exception='trypsin_exception'
        )
        canonical_peptides = set(canonical_peptides)

        with open(args.output_path, 'rt') as handle:
            novel_peptides = list(SeqIO.parse(handle, 'fasta'))

        # assert that the novel orf peptides don't overlap with the canonical
        # peptide
        self.assertEqual(
            len({str(x.seq) for x in novel_peptides}.intersection(canonical_peptides)),
            0
        )

        # assert that novel peptides are called from coding transcripts
        novel_orf_txs = set().union(
            *[{y.split('|')[0] for y in x.description.split(' ')} for x in novel_peptides]
        )
        # The check should verify that at least one transcript has coding transcripts
        # by checking if transcript IDs from novel peptides match any in proteome keys
        proteome_tx_ids = set(proteome.keys())
        coding_txs_found = novel_orf_txs.intersection(proteome_tx_ids)
        self.assertTrue(
            len(coding_txs_found) > 0,
            f"No coding transcripts found. Novel ORF txs: {novel_orf_txs}, Proteome txs: {proteome_tx_ids}"
        )
