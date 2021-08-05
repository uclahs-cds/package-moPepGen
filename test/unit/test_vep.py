""" Test the VEP data model """
import unittest
from test.unit import create_genomic_annotation, create_dna_record_dict
from moPepGen.parser import VEPParser


GENOME_DATA = {
    'chr1': 'ATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCT'
}
ANNOTATION_ATTRS = [
    {
        'gene_id': 'ENSG0001'
    },{
        'transcript_id': 'ENST0001.1',
        'gene_id': 'ENSG0001',
        'protein_id': 'ENSP0001'
    }
]
ANNOTTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 40, ANNOTATION_ATTRS[0]),
        'transcripts': ['ENST0001']
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

class TestVEPRecord(unittest.TestCase):
    """ Test case for VEP record """
    def test_vep_to_variant_record_case1(self):
        """ Test convert vep to variant record for SNV Tcc/Acc
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:10',
            allele='A',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='10',
            cds_position='10',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('Atg', 'Ctg'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'A')
        self.assertEqual(str(variant.alt), 'C')

    def test_vep_to_variant_record_case2(self):
        """ Test convert vep to variant record for SNV tCc/tTc
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTg', 'aCg'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 10)
        self.assertEqual(int(variant.location.end), 11)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'C')

    def test_vep_to_variant_record_case3(self):
        """ Test convert vep to variant record for INDEL tCc/tc
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('tAt', 'tt'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 11)
        self.assertEqual(str(variant.ref), 'TA')
        self.assertEqual(str(variant.alt), 'T')

    def test_vep_to_variant_record_case4(self):
        """ Test convert vep to variant record for INDEL tcc/tGcc cdna_position
        10. Pattern like this is saying that there is a insertion of G after
        the T at position 10.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='10',
            cds_position='10',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('tat', 'tGat'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TG')

    def test_vep_to_variant_record_case5(self):
        """ Test convert vep to variant record for INDEL TCC/- cdna_position
        10-12. Pattern like this is saying that there is a deletion of TCC
        from position 10 to 12.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='10-12',
            cds_position='10-12',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('TAT', '-'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 8)
        self.assertEqual(int(variant.location.end), 12)
        self.assertEqual(str(variant.ref), 'CTAT')
        self.assertEqual(str(variant.alt), 'C')

    def test_vep_to_variant_record_case6(self):
        """ Test convert vep to variant record for INDEL -/GAG cdna_position
        10-11. Pattern like this is saying that there is an insertion of GAG
        between position 10 to 11.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='10-11',
            cds_position='10-11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('-', 'GAG'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TGAG')

    def test_vep_to_variant_record_case7(self):
        """ Test convert vep to variant record for INDEL tcc/tGAGcc with
        cdna_position 10-11. Pattern like this is saying that there is an
        insertion of GAG between position 10 to 11.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='10-11',
            cds_position='10-11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('tat', 'tGAGat'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TGAG')

    def test_vep_to_variant_record_case8(self):
        """ Test convert vep to variant record with a deletion at the begining
        of the transcript sequence.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='1-3',
            cds_position='1-3',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('ctg', '-'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 0)
        self.assertEqual(int(variant.location.end), 3)
        self.assertEqual(str(variant.ref), 'CTG')
        self.assertEqual(str(variant.alt), 'TA')

    def test_vep_to_variant_record_case9(self):
        """ Test convert vep to variant record with a deletion at the begining
        of the transcript sequence but the transcript is cds_start_NF.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTTATION_DATA)
        anno.transcripts['ENST0001.1'].transcript.attributes['tag'] = ['cds_start_NF']
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:11',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='1-3',
            cds_position='1-3',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('ctg', '-'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(int(variant.location.start), 0)
        self.assertEqual(int(variant.location.end), 4)
        self.assertEqual(str(variant.ref), 'CTGG')
        self.assertEqual(str(variant.alt), 'G')


if __name__ == '__main__':
    unittest.main()
