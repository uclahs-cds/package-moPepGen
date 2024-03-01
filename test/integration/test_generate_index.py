""" Test the command line interface """
import argparse
import subprocess as sp
import sys
import copy
from unittest.mock import MagicMock
from test.integration import TestCaseIntegration
from moPepGen import cli, aa, params, err
from moPepGen.gtf import GenomicAnnotation, GenomicAnnotationOnDisk
from moPepGen.index import IndexDir, IndexMetadata
from moPepGen.version import MetaVersion


class TestGenerateIndex(TestCaseIntegration):
    """ Test cases for moPepGen generateIndex """
    def create_base_args(self) -> argparse.Namespace:
        """ create base args """
        args = argparse.Namespace()
        args.command = 'generateIndex'
        args.genome_fasta = self.data_dir / 'genome.fasta'
        args.annotation_gtf = self.data_dir / 'annotation.gtf'
        args.proteome_fasta = self.data_dir / 'translate.fasta'
        args.gtf_symlink = False
        args.reference_source = None
        args.invalid_protein_as_noncoding = False
        args.cleavage_rule = 'trypsin'
        args.cleavage_exception = 'trypsin_exception'
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.miscleavage = 2
        args.quiet = True
        args.output_dir = self.work_dir / 'index'
        return args

    def test_generate_index_cli(self):
        """ Test generateIndex cli """
        cmd = f"""
        {sys.executable} -m moPepGen.cli generateIndex \\
            -o {self.work_dir}/index \\
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

    def test_generate_index_case1(self):
        """ Test genreate index """
        args = self.create_base_args()
        cli.generate_index(args)
        files = {str(file.name) for file in args.output_dir.glob('*')}
        expected = {
            'genome.pkl', 'proteome.pkl',
            'annotation.gtf', 'annotation_gene.idx', 'annotation_tx.idx',
            'canonical_peptides.pkl', 'coding_transcripts.pkl',
            'metadata.json'
        }
        self.assertEqual(files, expected)

        index_dir = IndexDir(self.work_dir/'index')

        genome = index_dir.load_genome()
        self.assertTrue(len(genome) > 0)

        anno = index_dir.load_annotation()
        self.assertTrue(len(anno.genes) > 0)
        self.assertTrue(len(anno.transcripts) > 0)

        proteome = index_dir.load_proteome()
        self.assertTrue(len(proteome) > 0)

        canonical_peptides = index_dir.load_canonical_peptides()
        self.assertTrue(len(canonical_peptides) > 0)

        coding_tx = index_dir.load_coding_tx()
        self.assertTrue(len(coding_tx) > 0)


    def test_generate_index_case2(self):
        """ Test genreate index """
        args = self.create_base_args()
        args.gtf_symlink = True
        cli.generate_index(args)
        files = {str(file.name) for file in args.output_dir.glob('*')}
        expected = {
            'genome.pkl', 'proteome.pkl',
            'annotation.gtf', 'annotation_gene.idx', 'annotation_tx.idx',
            'canonical_peptides.pkl', 'coding_transcripts.pkl',
            'metadata.json'
        }
        self.assertEqual(files, expected)

class TestCaseGenomicAnnotationOnDisk(TestCaseIntegration):
    """ Test case for GenomicAnnotationOnDisk """
    def test_generate_index(self):
        """ Test generate index """
        proteome = aa.AminoAcidSeqDict()
        proteome.dump_fasta(self.data_dir/'translate.fasta')
        gtf_file = self.data_dir/'annotation.gtf'
        index_dir = IndexDir(self.work_dir)
        index_dir.save_annotation(gtf_file, proteome=proteome)
        expect = {
            'annotation.gtf', 'annotation_gene.idx', 'annotation_tx.idx'
        }
        received = {str(file.name) for file in self.work_dir.glob('*')}
        self.assertEqual(received, expect)

    def test_load_index(self):
        """ test load index """
        anno1 = GenomicAnnotation()
        anno1.dump_gtf(self.data_dir/'annotation.gtf')

        anno2 = GenomicAnnotationOnDisk()
        anno2.load_index(self.data_dir/'annotation.gtf', 'gencode')

        self.assertEqual(anno1.genes.keys(), anno2.genes.keys())
        self.assertEqual(anno1.transcripts.keys(), anno2.transcripts.keys())

class TestIndexDir(TestCaseIntegration):
    """ Test cases for IndexDir """
    def test_validate_metadata(self):
        """ Test validate metadata """
        index_dir = IndexDir(self.work_dir)
        index_version = MetaVersion()
        index_cleavage_params = params.CleavageParams(
            enzyme='Trypsin',
            exception='trypsin_exception',
            miscleavage=2,
            min_mw=500,
            min_length=7,
            max_length=25
        )
        metadata = IndexMetadata(
            version=index_version,
            cleavage_params=index_cleavage_params,
            source='test'
        )
        index_dir.load_metadata = MagicMock(return_value=metadata)
        index_dir.validate_metadata()
        self.assertTrue(index_dir.validate_metadata())

        cur_cleavage_params = copy.copy(index_cleavage_params)
        cur_cleavage_params.enzyme = 'LysC'
        cur_cleavage_params.exception = None
        with self.assertRaises(err.InvalidIndexError):
            index_dir.validate_metadata(cur_cleavage_params)
