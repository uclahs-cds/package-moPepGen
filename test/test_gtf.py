""" Test the GTF files are loaded and handled properly
"""
import unittest
from Bio import SeqIO
from moPepGen import gtf
from moPepGen.SeqFeature import SeqFeature, FeatureLocation


def create_transcript_model(data:dict) -> gtf.TranscriptAnnotationModel:
    """ Create a transcript model from data.
    """
    chrom = data['chrom']
    strand = data['strand']
    entry = data['transcript']
    location = FeatureLocation(start=entry[0], end=entry[1], strand=strand)
    transcript = SeqFeature(chrom=chrom, location=location,
        attributes=entry[2])
    exons = []
    for entry in data['exon']:
        location = FeatureLocation(start=entry[0], end=entry[1], strand=strand)
        exons.append(SeqFeature(chrom=chrom, location=location,
            attributes=entry[2]))
    cds = []
    if 'cds' in data:
        for entry in data['cds']:
            location = FeatureLocation(start=entry[0], end=entry[1],
                strand=strand)
            cds.append(SeqFeature(chrom=chrom, location=location,
                attributes=entry[2]))
    model = gtf.TranscriptAnnotationModel(transcript, cds, exons)
    return model

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
        chrom = SeqIO.read('test/files/genome_example.fa', 'fasta')
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
        chrom = SeqIO.read('test/files/genome_example.fa', 'fasta')
        seq = model.get_transcript_sequence(chrom)
        self.assertIs(seq.orf, None)

class TestGTF(unittest.TestCase):
    """ Test case for GTF modules
    """
    @staticmethod
    def load_gtf(path:str):
        """ Load the gtf file from disk """
        anno = gtf.TranscriptGTFDict()
        anno.dump_gtf(path)
        return anno

    def test_dump_gtf(self):
        """ Test that the VapIO.parse returns an iterabale of vep_record
        """
        anno = self.load_gtf('test/files/gtf_example_gencode.gtf')

        self.assertIsInstance(anno, gtf.TranscriptGTFDict)
        self.assertEqual(len(anno), 5)
        for key, val in anno.items():
            self.assertEqual(val.transcript.attributes['transcript_id'], key)
            for cds in val.cds:
                self.assertEqual(cds.attributes['transcript_id'], key)
            for exon in val.exon:
                self.assertEqual(exon.attributes['transcript_id'], key)


if __name__ == '__main__':
    unittest.main()
