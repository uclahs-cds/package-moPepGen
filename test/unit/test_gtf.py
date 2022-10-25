""" Test the GTF files are loaded and handled properly
"""
import copy
import io
import unittest
from test.unit import create_transcript_model, create_variant, \
    create_genomic_annotation
from Bio import SeqIO
from moPepGen import gtf, err
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from .test_vep_parser import ANNOTATION_DATA


ANNOTATION_ATTRS = [
    {
        'gene_id': 'ENSG0001'
    },{
        'transcript_id': 'ENST0001.1',
        'gene_id': 'ENSG0001',
        'protein_id': 'ENSP0001'
    }
]
ANNOTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 40, ANNOTATION_ATTRS[0]),
        'transcripts': ['ENST0001.1']
    }],
    'transcripts': [{
        # seq: CTGGT CCCCT ATGGG TCCTT C
        'transcript_id': ANNOTATION_ATTRS[1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (5, 35, ANNOTATION_ATTRS[1]),
        'exon': [
            (5, 12, ANNOTATION_ATTRS[1]),
            (17, 23, ANNOTATION_ATTRS[1]),
            (27, 35, ANNOTATION_ATTRS[1])
        ]
    }]
}

class TestGTFSeqFeature(unittest.TestCase):
    """ Test cases for GTFSeqFeature """
    def test_infer_source_ensembl(self):
        """ Test source is infered from ensembl """
        gtf_data = '1\thavana\ttranscript\t11869\t14409\t.\t+\t.\t' +\
            'gene_id "ENSG00000223972"; gene_version "5"; transcript_id'+\
            ' "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1";'+\
            ' gene_source "havana"; gene_biotype' +\
            ' "transcribed_unprocessed_pseudogene"; transcript_name'+\
            ' "DDX11L1-202"; transcript_source "havana"; transcript_biotype'+\
            ' "processed_transcript"; tag "basic"; transcript_support_level "1"'

        with io.BytesIO(gtf_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                record:GTFSeqFeature = list(gtf.GtfIO.parse(handle))[0]
        record.infer_annotation_source()
        self.assertEqual(record.source, 'ENSEMBL')
        self.assertEqual(record.biotype, 'transcribed_unprocessed_pseudogene')

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

    def test_get_transcript_sequence_case3(self):
        """ Test the transcript sequence with first CDS's frame isn't 0. """
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
        model.cds[0].frame = 1
        chrom = SeqIO.read('test/files/genome.fasta', 'fasta')
        seq = model.get_transcript_sequence(chrom)
        self.assertEqual(seq.orf.start, 1)

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
        self.assertEqual(i, 124)
        i = model.get_transcript_index(250)
        self.assertEqual(i, 74)

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
            self.assertEqual(val.transcript.transcript_id, key)
            for cds in val.cds:
                self.assertEqual(cds.transcript_id, key)
            for exon in val.exon:
                self.assertEqual(exon.transcript_id, key)

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
            location=FeatureLocation(seqname='chr22', start=0, end=350, strand=1),
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
            location=FeatureLocation(seqname='chr22', start=0, end=350, strand=1),
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

    def test_find_exon_index_case1(self):
        """ Test finding exon index with strand + """
        # pylint: disable=W0212
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx_id = 'ENST0001.1'
        exon = anno.transcripts[tx_id].exon[1]
        gene_id = 'ENSG0001'
        ind = anno.find_exon_index(tx_id, exon, 'genomic')
        self.assertEqual(ind, 1)

        exon2 = exon._shift(-anno.genes[gene_id].location.start)
        ind = anno.find_exon_index(tx_id, exon2, 'gene')
        self.assertEqual(ind, 1)

        exon = exon2._shift(10)
        with self.assertRaises(err.ExonNotFoundError):
            anno.find_exon_index(tx_id, exon, 'gene')

    def test_find_exon_index_case2(self):
        """ Test finding exon index with strand - """
        # pylint: disable=W0212

        anno = create_genomic_annotation(ANNOTATION_DATA)
        for gene in anno.genes.values():
            gene.location.strand = -1
        for transcript in anno.transcripts.values():
            transcript.transcript.location.strand = -1
            for exon in transcript.exon:
                exon.location.strand = -1

        tx_id = 'ENST0001.1'
        exon = anno.transcripts[tx_id].exon[0]
        gene_id = 'ENSG0001'
        ind = anno.find_exon_index(tx_id, exon, 'genomic')
        self.assertEqual(ind, 2)

        exon2 = anno.feature_coordinate_genomic_to_gene(exon, gene_id)
        ind = anno.find_exon_index(tx_id, exon2, 'gene')
        self.assertEqual(ind, 2)

        exon = exon2._shift(10)
        with self.assertRaises(err.ExonNotFoundError):
            anno.find_exon_index(tx_id, exon, 'gene')

    def test_find_intron_index_case1(self):
        """ Test finding intron index with strand = 1 """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001.1'
        start = anno.transcripts[tx_id].exon[0].location.end
        end = anno.transcripts[tx_id].exon[1].location.start
        strand = anno.genes[gene_id].location.strand
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=strand)
        intron = GTFSeqFeature(chrom=gene_id, location=location, attributes={})
        ind = anno.find_intron_index(tx_id, intron, 'genomic')
        self.assertEqual(ind, 0)

        intron2 = intron._shift(-anno.genes[gene_id].location.start)
        ind = anno.find_intron_index(tx_id, intron2, 'gene')
        self.assertEqual(ind, 0)

        intron = intron2._shift(10)
        with self.assertRaises(err.ExonNotFoundError):
            anno.find_exon_index(tx_id, intron, 'gene')

    def test_find_intron_index_case2(self):
        """ Test finding intron index with strand = -1 """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        for gene in anno.genes.values():
            gene.location.strand = -1
        for transcript in anno.transcripts.values():
            transcript.transcript.location.strand = -1
            for exon in transcript.exon:
                exon.location.strand = -1

        gene_id = 'ENSG0001'
        tx_id = 'ENST0001.1'
        start = anno.transcripts[tx_id].exon[0].location.end
        end = anno.transcripts[tx_id].exon[1].location.start
        strand = anno.genes[gene_id].location.strand
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=strand)
        intron = GTFSeqFeature(chrom=gene_id, location=location, attributes={})
        ind = anno.find_intron_index(tx_id, intron, 'genomic')
        self.assertEqual(ind, 1)

        intron2 = anno.feature_coordinate_genomic_to_gene(intron, gene_id)
        ind = anno.find_intron_index(tx_id, intron2, 'gene')
        self.assertEqual(ind, 1)

        intron = intron2._shift(10)
        with self.assertRaises(err.ExonNotFoundError):
            anno.find_exon_index(tx_id, intron, 'gene')

    def test_find_intron_index_fuzzy_1(self):
        """ Test finding intron index with end is earlier """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001.1'
        start = anno.transcripts[tx_id].exon[0].location.end
        end = anno.transcripts[tx_id].exon[1].location.start - 2
        strand = anno.genes[gene_id].location.strand
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=strand)
        intron = GTFSeqFeature(chrom=gene_id, location=location, attributes={})
        ind = anno.find_intron_index(tx_id, intron, 'genomic',
            intron_end_range=(-3, 0))
        self.assertEqual(ind, 0)

    def test_find_intron_index_fuzzy_2(self):
        """ Test finding intron index with star is earlier in - strand """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        for gene in anno.genes.values():
            gene.location.strand = -1
        for transcript in anno.transcripts.values():
            transcript.transcript.location.strand = -1
            for exon in transcript.exon:
                exon.location.strand = -1
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001.1'
        start = anno.transcripts[tx_id].exon[0].location.end
        end = anno.transcripts[tx_id].exon[1].location.start + 1
        strand = anno.genes[gene_id].location.strand
        location = FeatureLocation(seqname=gene_id, start=start, end=end,
            strand=strand)
        intron = GTFSeqFeature(chrom=gene_id, location=location, attributes={})
        ind = anno.find_intron_index(tx_id, intron, 'genomic',
            intron_start_range=(-2, 0))
        self.assertEqual(ind, 1)

    def test_coordinate_convert(self):
        """ Convert coodinates """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = ANNOTATION_ATTRS[0]['gene_id']
        tx_id = ANNOTATION_DATA['genes'][0]['transcripts'][0]
        x = anno.coordinate_gene_to_transcript(20, gene_id, tx_id)
        self.assertEqual(x, 10)

        annotation_data2 = copy.deepcopy(ANNOTATION_DATA)
        annotation_data2['genes'][0]['strand'] = -1
        annotation_data2['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(annotation_data2)
        gene_id = ANNOTATION_ATTRS[0]['gene_id']
        tx_id = ANNOTATION_DATA['genes'][0]['transcripts'][0]
        x = anno.coordinate_gene_to_transcript(19, gene_id, tx_id)
        self.assertEqual(x, 10)

    def test_coordinate_convert_tx_exon_start(self):
        """ Convert coodinate from transcript to genomic when the index is the
        start of an exon"""
        anno = create_genomic_annotation(ANNOTATION_DATA)
        gene_id = ANNOTATION_ATTRS[0]['gene_id']
        tx_id = ANNOTATION_DATA['genes'][0]['transcripts'][0]

        i_genomic = 17
        i_gene = anno.coordinate_genomic_to_gene(i_genomic, gene_id)
        i_tx = anno.coordinate_gene_to_transcript(i_gene, gene_id, tx_id)
        self.assertEqual(i_tx, 7)
        i_genomic_2 = anno.coordinate_transcript_to_genomic(i_tx, tx_id)
        self.assertEqual(i_genomic, i_genomic_2)

        i_genomic = 11
        i_gene = anno.coordinate_genomic_to_gene(i_genomic, gene_id)
        i_tx = anno.coordinate_gene_to_transcript(i_gene, gene_id, tx_id)
        self.assertEqual(i_tx, 6)
        i_genomic_2 = anno.coordinate_transcript_to_genomic(i_tx, tx_id)
        self.assertEqual(i_genomic, i_genomic_2)

        # negative strand
        annotation_data2 = copy.deepcopy(ANNOTATION_DATA)
        annotation_data2['genes'][0]['strand'] = -1
        annotation_data2['transcripts'][0]['strand'] = -1
        anno = create_genomic_annotation(annotation_data2)
        gene_id = ANNOTATION_ATTRS[0]['gene_id']
        tx_id = ANNOTATION_DATA['genes'][0]['transcripts'][0]

        i_genomic = 22
        i_gene = anno.coordinate_genomic_to_gene(i_genomic, gene_id)
        i_tx = anno.coordinate_gene_to_transcript(i_gene, gene_id, tx_id)
        self.assertEqual(i_tx, 8)
        i_genomic_2 = anno.coordinate_transcript_to_genomic(i_tx, tx_id)
        self.assertEqual(i_genomic, i_genomic_2)

        i_genomic = 17
        i_gene = anno.coordinate_genomic_to_gene(i_genomic, gene_id)
        i_tx = anno.coordinate_gene_to_transcript(i_gene, gene_id, tx_id)
        self.assertEqual(i_tx, 13)
        i_genomic_2 = anno.coordinate_transcript_to_genomic(i_tx, tx_id)
        self.assertEqual(i_genomic, i_genomic_2)
