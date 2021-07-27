""" Test the VEP data model """
import unittest
from Bio.Seq import Seq
from moPepGen import dna
from moPepGen.parser import VEPParser


class TestVEPRecord(unittest.TestCase):
    """ Test case for VEP record """
    def test_vep_to_variant_record_case1(self):
        """ Test convert vep to variant record for SNV Tcc/Acc
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
            codons=('Tcc', 'Acc'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'A')

    def test_vep_to_variant_record_case2(self):
        """ Test convert vep to variant record for SNV tCc/tTc
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
            codons=('tCc', 'tTc'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 10)
        self.assertEqual(int(variant.location.end), 11)
        self.assertEqual(str(variant.ref), 'C')
        self.assertEqual(str(variant.alt), 'T')

    def test_vep_to_variant_record_case3(self):
        """ Test convert vep to variant record for INDEL tCc/tc
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
            codons=('tCc', 'tc'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 11)
        self.assertEqual(str(variant.ref), 'TC')
        self.assertEqual(str(variant.alt), 'T')

    def test_vep_to_variant_record_case4(self):
        """ Test convert vep to variant record for INDEL tcc/tGcc cdna_position
        10. Pattern like this is saying that there is a insertion of G after
        the T at position 10.
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
            codons=('tcc', 'tGcc'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TG')

    def test_vep_to_variant_record_case5(self):
        """ Test convert vep to variant record for INDEL TCC/- cdna_position
        10-12. Pattern like this is saying that there is a deletion of TCC
        from position 10 to 12.
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
            codons=('TCC', '-'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 8)
        self.assertEqual(int(variant.location.end), 12)
        self.assertEqual(str(variant.ref), 'GTCC')
        self.assertEqual(str(variant.alt), 'G')

    def test_vep_to_variant_record_case6(self):
        """ Test convert vep to variant record for INDEL -/GAG cdna_position
        10-11. Pattern like this is saying that there is an insertion of GAG
        between position 10 to 11.
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TGAG')

    def test_vep_to_variant_record_case7(self):
        """ Test convert vep to variant record for INDEL tcc/tGAGcc with
        cdna_position 10-11. Pattern like this is saying that there is an
        insertion of GAG between position 10 to 11.
        """
        seq = dna.DNASeqRecord(
            seq=Seq('ATGTACTGGTCCTTCTGCCT')
        )
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
            codons=('tcc', 'tGAGcc'),
            existing_variation='-',
            extra={}
        )
        variant = vep_record.convert_to_variant_record(seq)
        self.assertEqual(int(variant.location.start), 9)
        self.assertEqual(int(variant.location.end), 10)
        self.assertEqual(str(variant.ref), 'T')
        self.assertEqual(str(variant.alt), 'TGAG')


if __name__ == '__main__':
    unittest.main()
