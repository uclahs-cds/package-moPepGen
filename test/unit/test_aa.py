""" Test module for amino acids """
import unittest
from Bio.Seq import Seq
from moPepGen import aa


class TestAminoAcidSeqDict(unittest.TestCase):
    """ Test case for AminoAcidSeqDict"""
    def testcreate_unique_peptide_pool_by_length(self):
        """ Test that peptides kept are of correct length. """
        sequence = 'MKKVTAEAISWNESTSETNNSMVTEFIFLGLSDSQELQTFLFMLFFVFYGGIVFGN'+\
            'LLIVITVVSDSHLHSPMYFLLANLSLIDLSLSSVTAPKMITDFFSQRKVISFKGCLVQIFLLH'+\
            'FFGGSEMVILIAMGFDRYIAICKPLHYTTIMCGNACVGIMAVTWGIGFLHSVSQLAFAVHLLF'+\
            'CGPNEVDSFYCDLPRVIKLACTDTYRLDIMVIANSGVLTVCSFVLLIISYTIILMTIQHRPLD'+\
            'KSSKALSTLTAHITVVLLFFGPCVFIYAWPFPIKSLDKFLAVFYSVITPLLNPIIYTLRNKDM'+\
            'KTAIRQLRKWDAHSSVKF'
        header = 'ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.6|OTTH'+\
            'UMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|326'
        seq = aa.AminoAcidSeqRecord(
            seq=sequence,
            _id=header,
            name=header,
            description=header
        )
        seq.infer_ids()
        kwargs = {seq.transcript_id : seq}
        seqdict = aa.AminoAcidSeqDict(**kwargs)
        pool = seqdict.create_unique_peptide_pool(rule='trypsin',
            exception='trypsin_exception', miscleavage=2, min_mw=0.,
            min_length = 7, max_length = 25)
        expected = {'TAIRQLR', 'WDAHSSVKF', 'MITDFFSQR',
            'SLDKFLAVFYSVITPLLNPIIYTLR', 'NKDMKTAIR', 'MITDFFSQRK',
            'TAIRQLRK', 'LACTDTYR', 'DMKTAIRQLR', 'QLRKWDAHSSVK', 'DMKTAIR',
            'FLAVFYSVITPLLNPIIYTLRNK', 'MITDFFSQRKVISFK', 'WDAHSSVK',
            'KWDAHSSVKF', 'VIKLACTDTYR', 'KWDAHSSVK', 'FLAVFYSVITPLLNPIIYTLR'}
        itol = {'TALRQLR', 'MLTDFFSQR', 'SLDKFLAVFYSVLTPLLNPLLYTLR',
            'NKDMKTALR', 'MLTDFFSQRK', 'TALRQLRK', 'DMKTALRQLR', 'DMKTALR',
            'FLAVFYSVLTPLLNPLLYTLRNK', 'MLTDFFSQRKVLSFK',
            'VLKLACTDTYR', 'FLAVFYSVLTPLLNPLLYTLR'}
        expected_pool = expected.union(itol)
        self.assertEqual(pool, expected_pool)

class TestAminoAcidSeqRecord(unittest.TestCase):
    """ Test case for AminoAcidSeqRecord """
    def testinfer_ids_ensembl_case1(self):
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
        self.assertEqual(seq.protein_id, 'ENSP00000488240.1')
        self.assertEqual(seq.transcript_id, 'ENST00000631435.1')
        self.assertEqual(seq.gene_id, 'ENSG00000282253.1')

    def testinfer_ids_ensembl_case2(self):
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

    def testinfer_ids_gencode_case1(self):
        """ Test that ids are infered correctly with ENSMBLE style. """
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

    def testinfer_ids_gencode_case2(self):
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
