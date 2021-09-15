""" Test the VEP data model """
import copy
import unittest
from test.unit import create_genomic_annotation, create_dna_record_dict
from moPepGen.parser import VEPParser
from moPepGen.err import TranscriptionStopSiteMutationError, \
    TranscriptionStartSiteMutationError


GENOME_DATA = {
    'chr1': 'ATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCT'
}
ANNOTATION_ATTRS = [
    {
        'gene_id': 'ENSG0001',
        'gene_name': 'SYMBO'
    },{
        'transcript_id': 'ENST0001.1',
        'gene_id': 'ENSG0001',
        'protein_id': 'ENSP0001',
        'gene_name': 'SYMBO'
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

class TestVEPRecord(unittest.TestCase):
    """ Test case for VEP record """
    def test_vep_to_variant_record_case1(self):
        """ Test convert vep to variant record for SNV Tcc/Acc
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:21',
            allele='C',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('Atg', 'Ctg'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        tx_variant = variant.to_transcript_variant(anno, genome)
        self.assertEqual(int(tx_variant.location.start), 10)
        self.assertEqual(int(tx_variant.location.end), 11)
        self.assertEqual(str(tx_variant.ref), 'A')
        self.assertEqual(str(tx_variant.alt), 'C')

    def test_vep_to_variant_record_case2(self):
        """ Test convert vep to variant record for SNV tCc/tTc
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:22',
            allele='C',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='12',
            cds_position='12',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTg', 'aCg'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        variant = variant.to_transcript_variant(anno, genome)
        self.assertEqual(int(variant.location.start), 11)
        self.assertEqual(int(variant.location.end), 12)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'C')

    def test_vep_to_variant_record_case3(self):
        """ Test convert vep to variant record for INDEL tCc/tc
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:21',
            allele='-',
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
        variant = variant.to_transcript_variant(anno, genome)
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
        anno = create_genomic_annotation(ANNOTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:20-21',
            allele='G',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='10-11',
            cds_position='10-11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('tat', 'tGat'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        variant = variant.to_transcript_variant(anno, genome)
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
        anno = create_genomic_annotation(ANNOTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:20-22',
            allele='-',
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
        variant = variant.to_transcript_variant(anno, genome)
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
        anno = create_genomic_annotation(ANNOTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:20-21',
            allele='GAG',
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
        variant = variant.to_transcript_variant(anno, genome)
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
        anno = create_genomic_annotation(ANNOTATION_DATA)
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:20-21',
            allele='GAG',
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
        variant = variant.to_transcript_variant(anno, genome)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TGAG')

    def test_vep_to_variant_record_case9(self):
        """ Test convert vep to variant record with a deletion at the begining
        of the transcript sequence but the transcript is cds_start_NF.
        """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)
        anno.transcripts['ENST0001.1'].transcript.attributes['tag'] = ['cds_start_NF']
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:6-8',
            allele='-',
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
        variant = variant.to_transcript_variant(anno, genome)
        self.assertEqual(int(variant.location.start), 0)
        self.assertEqual(int(variant.location.end), 4)
        self.assertEqual(str(variant.ref), 'CTGG')
        self.assertEqual(str(variant.alt), 'G')

    def test_vep_to_variant_record_case10(self):
        """ Test convert vep to variant record with - strand
        """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:21',
            allele='C',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTa', 'aCa'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(anno, genome)
        tx_variant = variant.to_transcript_variant(anno, genome)
        self.assertEqual(int(tx_variant.location.start), 10)
        self.assertEqual(int(tx_variant.location.end), 11)
        self.assertEqual(str(tx_variant.ref), 'T')
        self.assertEqual(str(tx_variant.alt), 'G')

    def test_vep_to_variant_record_case11(self):
        """ Transcriptional stop altering variant """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:38-50',
            allele='-',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTa', 'aCa'),
            existing_variation='-',
            extra={}
        )
        with self.assertRaises(TranscriptionStopSiteMutationError):
            vep_record.convert_to_variant_record(anno, genome)

    def test_vep_to_variant_record_case12(self):
        """ Transcriptional stop altering variant on - strand """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:3-8',
            allele='-',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTa', 'aCa'),
            existing_variation='-',
            extra={}
        )
        with self.assertRaises(TranscriptionStopSiteMutationError):
            vep_record.convert_to_variant_record(anno, genome)

    def test_vep_to_variant_record_case13(self):
        """ Transcriptional stop altering variant """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:10-15',
            allele='-',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTa', 'aCa'),
            existing_variation='-',
            extra={}
        )
        record = vep_record.convert_to_variant_record(anno, genome)
        res = record.is_spanning_over_splicing_site(anno, 'ENST0001.1')
        self.assertTrue(res)

    def test_vep_to_variant_record_case14(self):
        """ Test convert vep to variant record with a deletion at the begining
        of the transcript sequence on - strand.
        """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno_data['genes'][0]['strand'] = -1
        anno_data['transcripts'][0]['strand'] = -1
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(anno_data)
        anno.transcripts['ENST0001.1'].transcript.attributes['tag'] = ['cds_start_NF']
        # seq: CTGGT CCCCT ATGGG TCCTT C
        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:33-35',
            allele='-',
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
        variant = variant.to_transcript_variant(anno, genome)
        self.assertEqual(int(variant.location.start), 0)
        self.assertEqual(int(variant.location.end), 4)
        self.assertEqual(str(variant.ref), 'GAAG')
        self.assertEqual(str(variant.alt), 'G')

    def test_vep_to_variant_record_case15(self):
        """ Transcriptional start altering variant """
        annotation_data = copy.copy(ANNOTATION_DATA)
        annotation_data['genes'][0]['gene'] = (5, 40, ANNOTATION_ATTRS[0])
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(annotation_data)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:2-7',
            allele='-',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTa', 'aCa'),
            existing_variation='-',
            extra={}
        )
        with self.assertRaises(TranscriptionStartSiteMutationError):
            vep_record.convert_to_variant_record(anno, genome)

    def test_vep_to_variant_record_case15_deletion(self):
        """ deletion that allele is not - """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        vep_record = VEPParser.VEPRecord(
            uploaded_variation='rs55971985',
            location='chr1:19-22',
            allele='T',
            gene='ENSG0001',
            feature='ENST0001.1',
            feature_type='Transcript',
            consequences=['missense_variant'],
            cdna_position='11',
            cds_position='11',
            protein_position=3,
            amino_acids=('S', 'T'),
            codons=('aTa', 'aCa'),
            existing_variation='-',
            extra={}
        )
        record = vep_record.convert_to_variant_record(anno, genome)
        self.assertEqual(record.ref, 'CTAT')
        self.assertEqual(record.alt, 'T')

if __name__ == '__main__':
    unittest.main()
