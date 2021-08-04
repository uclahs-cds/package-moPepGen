""" Test the GTF files are loaded and handled properly
"""
import unittest
from test.unit import create_transcript_model, create_variant
from Bio import SeqIO
from moPepGen import gtf
from moPepGen.SeqFeature import FeatureLocation


class TestAnnotationModel(unittest.TestCase):
    """ Test case for the annotation model """
    def test_get_transcript_sequence_case1(self):
        """ Test the transcript sequence is returned correctly. """
        attributes = {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001'
        }
        data = {
            'chrom': 'chr22',
            'strand': 1,
            'transcript': (0, 20, attributes),
            'exon': [(0,8,attributes), (10, 18, attributes)],
            'cds': [(0,8,attributes), (10, 18, attributes)]
        }
        model = create_transcript_model(data)
        chrom = SeqIO.read('test/files/genome.fasta', 'fasta')
        seq = model.get_transcript_sequence(chrom)
        self.assertEqual(len(seq.seq), 16)

    def test_get_transcript_sequence_case2(self):
        """ When there is no cds, the orf of the returned sequence should be
        None """
        attributes = {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001'
        }
        data = {
            'chrom': 'chr22',
            'strand': 1,
            'transcript': (0, 20, attributes),
            'exon': [(0,8,attributes), (10, 18, attributes)]
        }
        model = create_transcript_model(data)
        chrom = SeqIO.read('test/files/genome.fasta', 'fasta')
        seq = model.get_transcript_sequence(chrom)
        self.assertIs(seq.orf, None)

    def test_get_transcript_index_case1(self):
        """ Getting transcript index from genomic index when strand is + """
        attributes = {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001'
        }
        data = {
            'chrom': 'chr22',
            'strand': 1,
            'transcript': (150, 350, attributes),
            'exon': [
                (150,200,attributes),
                (225, 275, attributes),
                (300, 350, attributes)
            ]
        }
        model = create_transcript_model(data)
        i = model.get_transcript_index(175)
        self.assertEqual(i, 25)
        i = model.get_transcript_index(250)
        self.assertEqual(i, 75)

        with self.assertRaises(ValueError):
            model.get_transcript_index(0)

        with self.assertRaises(ValueError):
            model.get_transcript_index(210)

        with self.assertRaises(ValueError):
            model.get_transcript_index(375)

    def test_get_transcript_index_case2(self):
        """ Getting transcript index from genomic index when strand is - """
        attributes = {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001'
        }
        data = {
            'chrom': 'chr22',
            'strand': -1,
            'transcript': (150, 350, attributes),
            'exon': [
                (150,200,attributes),
                (225, 275, attributes),
                (300, 350, attributes)
            ]
        }
        model = create_transcript_model(data)
        i = model.get_transcript_index(175)
        self.assertEqual(i, 125)
        i = model.get_transcript_index(250)
        self.assertEqual(i, 75)

        with self.assertRaises(ValueError):
            model.get_transcript_index(0)

        with self.assertRaises(ValueError):
            model.get_transcript_index(210)

        with self.assertRaises(ValueError):
            model.get_transcript_index(375)

class TestGTF(unittest.TestCase):
    """ Test case for GTF modules
    """
    @staticmethod
    def load_gtf(path:str):
        """ Load the gtf file from disk """
        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(path)
        return anno

    def test_dump_gtf(self):
        """ Test that the VapIO.parse returns an iterabale of vep_record
        """
        anno = self.load_gtf('test/files/annotation.gtf')

        self.assertIsInstance(anno, gtf.GenomicAnnotation)
        self.assertEqual(len(anno.transcripts), 5)
        for key, val in anno.transcripts.items():
            self.assertEqual(val.transcript.attributes['transcript_id'], key)
            for cds in val.cds:
                self.assertEqual(cds.attributes['transcript_id'], key)
            for exon in val.exon:
                self.assertEqual(exon.attributes['transcript_id'], key)

    def test_variant_coordinates_convert_case1(self):
        """ Test the converting the coordinates of variants to gene """
        attributes = {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001'
        }
        data = {
            'chrom': 'chr22',
            'strand': 1,
            'transcript': (150, 350, attributes),
            'exon': [
                (150,200,attributes),
                (225, 275, attributes),
                (300, 350, attributes)
            ]
        }
        transcript_model = create_transcript_model(data)
        transcript_model.transcript.strand = 1
        gene_model = gtf.GeneAnnotationModel(
            location=FeatureLocation(seqname='chr22', start=0, end=350),
            chrom='chr22',
            transcripts=[attributes['transcript_id']],
            attributes={}
        )
        genes = {attributes['gene_id']: gene_model}
        transcripts = {attributes['transcript_id']: transcript_model}
        anno = gtf.GenomicAnnotation(genes, transcripts)
        var1 = create_variant(25, 26, 'A', 'C', 'SNV', 'XX')
        var1.location.seqname = attributes['transcript_id']
        var2 = anno.variant_coordinates_to_gene(var1, attributes['gene_id'])
        self.assertEqual(int(var2.location.start), 175)
        self.assertEqual(int(var2.location.end), 176)

    def test_variant_coordinates_convert_case2(self):
        """ Test the converting the coordinates of variants to gene """
        attributes = {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001'
        }
        data = {
            'chrom': 'chr22',
            'strand': -1,
            'transcript': (150, 350, attributes),
            'exon': [
                (150,200,attributes),
                (225, 275, attributes),
                (300, 350, attributes)
            ]
        }
        transcript_model = create_transcript_model(data)
        transcript_model.transcript.strand = 1
        gene_model = gtf.GeneAnnotationModel(
            location=FeatureLocation(seqname='chr22', start=0, end=350),
            chrom='chr22',
            transcripts=[attributes['transcript_id']],
            attributes={}
        )
        genes = {attributes['gene_id']: gene_model}
        transcripts = {attributes['transcript_id']: transcript_model}
        anno = gtf.GenomicAnnotation(genes, transcripts)
        var1 = create_variant(25, 26, 'A', 'C', 'SNV', 'XX')
        var1.location.seqname = attributes['transcript_id']
        var2 = anno.variant_coordinates_to_gene(var1, attributes['gene_id'])
        self.assertEqual(int(var2.location.start), 175)
        self.assertEqual(int(var2.location.end), 176)


if __name__ == '__main__':
    unittest.main()