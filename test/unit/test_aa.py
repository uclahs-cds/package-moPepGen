""" Test module for amino acids """
import unittest
from test.unit import create_aa_record, create_genomic_annotation, \
    create_aa_seq_with_coordinates
from Bio.Seq import Seq
from moPepGen import aa, params


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
        'transcripts': ['ENST0001']
    }],
    'transcripts': [{
        'transcript_id': ANNOTATION_ATTRS[1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (5, 35, ANNOTATION_ATTRS[1]),
        'exon': []
    }]
}

class TestAminoAcidSeqDict(unittest.TestCase):
    """ Test case for AminoAcidSeqDict"""
    def test_create_unique_peptide_pool_by_length(self):
        """ Test that peptides kept are of correct length. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        sequence = 'MKVTAEAISWNERSTSETNNSMVTEFIFLGLSDSQEKLQRTFLFMLFFVFYGGIVFGN'
        header = 'ENSP0001|ENST0001.1|ENSG0001|-'
        seq = aa.AminoAcidSeqRecord(
            seq=sequence,
            _id=header,
            name=header,
            description=header
        )
        seq.infer_ids()
        kwargs = {seq.transcript_id : seq}
        seqdict = aa.AminoAcidSeqDict(**kwargs)
        pool = seqdict.create_unique_peptide_pool(anno=anno, rule='trypsin',
            exception='trypsin_exception', miscleavage=2, min_mw=0.,
            min_length = 7, max_length = 25)
        expected = {'MKVTAEAISWNER', 'KVTAEAISWNER', 'VTAEAISWNER',
            'STSETNNSMVTEFIFLGLSDSQEK', 'LQRTFLFMLFFVFYGGIVFGN',
            'TFLFMLFFVFYGGIVFGN', 'MKVTAEALSWNER', 'KVTAEALSWNER',
            'VTAEALSWNER', 'STSETNNSMVTEFLFLGLSDSQEK', 'LQRTFLFMLFFVFYGGLVFGN',
            'TFLFMLFFVFYGGLVFGN'}
        self.assertEqual(pool, expected)

class TestAminoAcidSeqRecord(unittest.TestCase):
    """ Test case for AminoAcidSeqRecord """
    def test_infer_ids_ensembl_case1(self):
        """ Test that ids are infered correctly with ENSMBLE style. """
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id='ENSP00000488240.1',
            name='ENSP00000488240.1',
            description='ENSP00000488240.1 pep chromosome:GRCh38:CHR_HSCHR7_2_'
            'CTG6:142847306:142847317:1 gene:ENSG00000282253.1 transcript:ENST'
            '00000631435.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene'
            ' gene_symbol:TRBD1 description:T cell receptor beta diversity 1 ['
            'Source:HGNC Symbol;Acc:HGNC:12158]'
        )
        seq.infer_ids_ensembl()
        self.assertEqual(seq.protein_id, 'ENSP00000488240')
        self.assertEqual(seq.transcript_id, 'ENST00000631435')
        self.assertEqual(seq.gene_id, 'ENSG00000282253')

    def test_infer_ids_ensembl_case2(self):
        """ Test that error will raise with GENCODE style. """
        header = 'ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.6|OTTH'+\
            'UMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|326'
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id=header,
            name=header,
            description=header
        )
        with self.assertRaises(ValueError):
            seq.infer_ids_ensembl()

    def test_infer_ids_ensembl_case3(self):
        """ Test that IDs are infered correctly with ENSEMBL style for GRCh37. """
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id='ENSP00000488240.1',
            name='ENSP00000488240.1',
            description='ENSP00000488240 pep:known chromosome:GRCh37:CHR_HSCHR7_2_'
            'CTG6:142847306:142847317:1 gene:ENSG00000282253 transcript:ENST'
            '00000631435 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene'
            ' gene_symbol:TRBD1 description:T cell receptor beta diversity 1 ['
            'Source:HGNC Symbol;Acc:HGNC:12158]'
        )
        seq.infer_ids_ensembl()
        self.assertEqual(seq.protein_id, 'ENSP00000488240')
        self.assertEqual(seq.transcript_id, 'ENST00000631435')
        self.assertEqual(seq.gene_id, 'ENSG00000282253')

    def test_infer_ids_gencode_case1(self):
        """ Test that ids are infered correctly with GENCODE style. """
        header = 'ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.6|OTTH'+\
            'UMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|326'
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id=header,
            name=header,
            description=header
        )
        seq.infer_ids_gencode()
        self.assertEqual(seq.id, 'ENSP00000493376.2')
        self.assertEqual(seq.name, 'ENSP00000493376.2')
        self.assertEqual(seq.protein_id, 'ENSP00000493376.2')
        self.assertEqual(seq.gene_id, 'ENSG00000186092.6')
        self.assertEqual(seq.transcript_id, 'ENST00000641515.2')

    def test_infer_ids_gencode_case2(self):
        """ Test that error will raise with ENSEMBL style """
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id='ENSP00000488240.1',
            name='ENSP00000488240.1',
            description='ENSP00000488240.1 pep chromosome:GRCh38:CHR_HSCHR7_2_'
            'CTG6:142847306:142847317:1 gene:ENSG00000282253.1 transcript:ENST'
            '00000631435.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene'
            ' gene_symbol:TRBD1 description:T cell receptor beta diversity 1 ['
            'Source:HGNC Symbol;Acc:HGNC:12158]'
        )
        with self.assertRaises(ValueError):
            seq.infer_ids_gencode()

    def test_infer_ids_gencode_case3(self):
        """ Test that ids are infered correctly with GENCODE style with _PAR_Y. """
        header = 'ENSP00000493376.2|ENST00000641515.2_PAR_Y|ENSG00000186092.6_PAR_Y' +\
            '|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|326'
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id=header,
            name=header,
            description=header
        )
        seq.infer_ids_gencode()
        self.assertEqual(seq.id, 'ENSP00000493376.2')
        self.assertEqual(seq.name, 'ENSP00000493376.2')
        self.assertEqual(seq.protein_id, 'ENSP00000493376.2')
        self.assertEqual(seq.gene_id, 'ENSG00000186092.6_PAR_Y')
        self.assertEqual(seq.transcript_id, 'ENST00000641515.2_PAR_Y')

    def test_auto_infer_ids_gencode_case1(self):
        """ Test that ids are infered correctly with GENCODE style """
        header = 'ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.6|OTTH'+\
            'UMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|326'
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id=header,
            name=header,
            description=header
        )
        source = seq.infer_ids()
        self.assertEqual(source, 'GENCODE')

    def test_auto_infer_ids_ensembl_case1(self):
        """ Test that ids are infered correctly with ENSMBLE style. """
        seq = aa.AminoAcidSeqRecord(
            seq=Seq('GTGG'),
            _id='ENSP00000488240.1',
            name='ENSP00000488240.1',
            description='ENSP00000488240.1 pep chromosome:GRCh38:CHR_HSCHR7_2_'
            'CTG6:142847306:142847317:1 gene:ENSG00000282253.1 transcript:ENST'
            '00000631435.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene'
            ' gene_symbol:TRBD1 description:T cell receptor beta diversity 1 ['
            'Source:HGNC Symbol;Acc:HGNC:12158]'
        )
        source = seq.infer_ids()
        self.assertEqual(source, 'ENSEMBL')

    def test_enzyme_lysn(self):
        """ Ensures that lysN cleaves lysine at N-terminus """
        seq = aa.AminoAcidSeqRecord('ACDEGKILMNP')
        expected = {'ACDEG', 'KILMNP'}
        fragments = seq.enzymatic_cleave(rule='lysn', miscleavage=0, min_mw=0, min_length=0)
        received = {str(x.seq) for x in fragments}
        self.assertEqual(expected, received)

class TestCaseVariantPeptidePool(unittest.TestCase):
    """ Test cases for VariantPeptidePool """
    def test_variant_peptide_add(self):
        """ Test variant peptide is added to the pool """
        peptide_data = [
            ('MFAEHTPK', 'ENSG0001|SNV-100-T-C|1'),
            ('GLAKER', 'ENSG0002|SNV-100-T-C|1'),
            ('MFAEHTPK', 'ENSG0003|SNV-100-T-C|1')
        ]
        peptides = [create_aa_record(*x) for x in peptide_data]
        pool = aa.VariantPeptidePool(set(peptides[:2]))
        canonical = {'ABCD'}
        # pool2 = {1,2,3}
        # get_equivalent(pool2, 2)p
        cleavage_params = params.CleavageParams()
        pool.add_peptide(peptides[-1], canonical, cleavage_params)
        for seq in pool.peptides:
            if str(seq.seq) == peptide_data[-1][0]:
                expected = f'{peptide_data[0][1]} {peptide_data[-1][1]}'
                self.assertEqual(seq.description, expected)

class TestCaseAminoAcidSeqRecordWithCoordinates(unittest.TestCase):
    """ Test cases for AminoAcidSeqRecordWithCoordinates """
    def test_add(self):
        """ test add """
        loc1 = [((0,4),(0,4))]
        loc2 = [((0,4),(4,8))]
        seq1 = create_aa_seq_with_coordinates('SSSR', loc1, (0,None))
        seq2 = create_aa_seq_with_coordinates('SSSF', loc2, (0,None))
        seq = seq1 + seq2
        self.assertEqual(str(seq.seq), 'SSSRSSSF')
        received = {((int(x.query.start), int(x.query.end)), \
            (int(x.ref.start), int(x.ref.end))) for x in seq.locations}
        expected = {((0,8),(0,8))}
        self.assertEqual(received, expected)

    def test_getitem(self):
        """ test getitem """
        loc1 = [((0,4),(0,4))]
        seq1 = create_aa_seq_with_coordinates('SSSR', loc1, (0,None))
        seq = seq1[1:3]
        self.assertEqual(str(seq.seq), 'SS')
        self.assertEqual(len(seq.locations), 1)
        self.assertEqual(seq.locations[0].query.start, 0)
        self.assertEqual(seq.locations[0].query.end, 2)
        self.assertEqual(seq.locations[0].ref.start, 1)
        self.assertEqual(seq.locations[0].ref.end, 3)

        loc1 = [((0,4),(0,4)), ((5,8), (5,8))]
        seq1 = create_aa_seq_with_coordinates('SSSRSSSG', loc1, (0,None))
        seq = seq1[3:7]
        self.assertEqual(len(seq.locations), 2)
        self.assertEqual(seq.locations[0].query.start, 0)
        self.assertEqual(seq.locations[0].query.end, 1)
        self.assertEqual(seq.locations[0].ref.start, 3)
        self.assertEqual(seq.locations[0].ref.end, 4)
        self.assertEqual(seq.locations[1].query.start, 2)
        self.assertEqual(seq.locations[1].query.end, 4)
        self.assertEqual(seq.locations[1].ref.start, 5)
        self.assertEqual(seq.locations[1].ref.end, 7)

        seq = seq1[4:7]
        self.assertEqual(len(seq.locations), 1)
        self.assertEqual(seq.locations[0].query.start, 1)
        self.assertEqual(seq.locations[0].query.end, 3)
        self.assertEqual(seq.locations[0].ref.start, 5)
        self.assertEqual(seq.locations[0].ref.end, 7)

        seq = seq1[2:5]
        self.assertEqual(len(seq.locations), 1)
        self.assertEqual(seq.locations[0].query.start, 0)
        self.assertEqual(seq.locations[0].query.end, 2)
        self.assertEqual(seq.locations[0].ref.start, 2)
        self.assertEqual(seq.locations[0].ref.end, 4)
