""" Module for testing PeptidePoolSplitter """
import copy
import io
import unittest
from test.unit import create_aa_record
from moPepGen.aa import PeptidePoolSplitter, VariantSourceSet, VariantPeptideInfo


GVF_CASE1 = [
    ['##source=gSNP'],
    ['##chrom=GENE_ID'],
    ['##parser=parseVEP'],
    ['ENSG0001', '611', 'SNV-611-G-A','G', 'A', '.', '.', 'TRANSCRIPT_ID=ENST0001']
]

PEPTIDE_DB_CASE1 = [
    '>ENST0001|SNV-1157-G-A|6',
    'IQGARATLLFATLLFER',
    '>ENST0002|SNV-100-G-T|3 ENST0003|INDEL-200-TGAAG-T|1',
    'MEGAGPRGAVPARR',
    '>ENST0004|SNV-100-G-T|7',
    'GAVPARRR'
]

LABEL_MAP1 = {
    'ENST0001': { 'SNV-1157-G-A': 'gSNP' },
    'ENST0002': { 'SNV-100-G-T': 'gSNP' },
    'ENST0003': { 'INDEL-200-TGAAG-T': 'gINDEL' },
    'ENST0004': { 'SNV-100-G-T': 'gSNP' }
}

class TestVariantSourceSet(unittest.TestCase):
    """ Test cases for VariantSourceSet """

    def setUp(self):
        """ Set up """
        super().setUp()
        VariantSourceSet.reset_levels()

    def test_source_levels(self):
        """ Test comparing variant source set """
        levels = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3}
        VariantSourceSet.set_levels(levels)
        set1 = VariantSourceSet(['gSNP'])
        self.assertEqual(set1.levels, levels)

        levels2 = copy.copy(levels)
        levels2['Fusion'] = 4
        set1.set_levels(levels2)
        set2 = VariantSourceSet(['sSNV'])
        self.assertEqual(set2.levels, levels2)

    def test_comparison(self):
        """ Test comparisons """
        levels = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3}
        VariantSourceSet.set_levels(levels)
        set1 = VariantSourceSet(['gSNP'])
        set2 = VariantSourceSet(['sSNV'])
        self.assertTrue(set1 < set2)

        set2 = VariantSourceSet(['gSNP', 'gINDEL'])
        self.assertTrue(set1 < set2)

        set1 = VariantSourceSet(['gINDEL'])
        self.assertTrue(set1 < set2)

        set1 = VariantSourceSet(['gINDEL', 'sSNV'])
        set2 = VariantSourceSet(['gINDEL', 'sINDEL'])
        self.assertTrue(set1 < set2)

class TestVariantPeptideInfo(unittest.TestCase):
    """ Test VariantPeptideInfo """
    def test_from_variant_peptide(self):
        """ Test forming from variant peptide """
        levels = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3}
        VariantSourceSet.set_levels(levels)
        label_map = {'ENST0001': {'SNV-1157-G-A': 'sSNV'}}
        peptide = create_aa_record('KHIRJ','ENST0001|SNV-1157-G-A|1')
        infos = VariantPeptideInfo.from_variant_peptide(peptide, label_map)
        self.assertEqual([x.sources for x in infos], [{'sSNV'}])

        label_map = {
            'ENST0001': {'SNV-1157-G-A': 'gSNP', 'INDEL-100-GGGGG-G': 'gINDEL'}
        }
        label = 'ENST0001|SNV-1157-G-A|INDEL-100-GGGGG-G|1'
        peptide = create_aa_record('KHIRJ', label)
        infos = VariantPeptideInfo.from_variant_peptide(peptide, label_map)
        self.assertEqual([x.sources for x in infos], [{'gSNP', 'gINDEL'}])

        label_map = {
            'ENST0001': {'SNV-1157-G-A': 'gSNP', 'INDEL-100-GGGGG-G': 'gINDEL'},
            'ENST0002': {'SNV-254-T-A': 'sSNV'}
        }
        label = 'ENST0001|SNV-1157-G-A|INDEL-100-GGGGG-G|1 ENST0002|SNV-254-T-A|1'
        peptide = create_aa_record('KHIRJ', label)
        infos = VariantPeptideInfo.from_variant_peptide(peptide, label_map)
        self.assertEqual([x.sources for x in infos], [{'gSNP', 'gINDEL'}, {'sSNV'}])

class TestPeptidePoolSplitter(unittest.TestCase):
    """ Test cases for testing PeptidePoolSplitter """
    def test_append_order_noncoding(self):
        """ Test the group order noncoding """
        levels = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3}
        splitter = PeptidePoolSplitter(order=levels)
        splitter.append_order_noncoding()
        self.assertEqual(splitter.order['Noncoding'], 4)

    def test_load_gvf(self):
        """ test loading gvf """
        levels = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3}
        splitter = PeptidePoolSplitter(order=levels)

        file_data = '\n'.join(['\t'.join(x) for x in GVF_CASE1])
        with io.BytesIO(file_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_gvf(handle)
        self.assertEqual(
            splitter.label_map['ENST0001']['SNV-611-G-A'], 'gSNP'
        )

    def test_load_gvf_with_unspecified_group(self):
        """ test loading gvf that the group is not specified. """
        levels = {'gINDEL': 0}
        splitter = PeptidePoolSplitter(order=levels)

        file_data = '\n'.join(['\t'.join(x) for x in GVF_CASE1])
        with io.BytesIO(file_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_gvf(handle)
        self.assertEqual(
            splitter.label_map['ENST0001']['SNV-611-G-A'], 'gSNP'
        )
        self.assertEqual(splitter.order['gSNP'], 1)

    def test_load_database(self):
        """ Test loading peptide database """
        splitter = PeptidePoolSplitter()

        db_data = '\n'.join(PEPTIDE_DB_CASE1)

        with io.BytesIO(db_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_database(handle)

        self.assertTrue(len(splitter.peptides.peptides) > 0)

    def test_load_database_with_noncoding(self):
        """ Test loading database of noncoding """
        splitter = PeptidePoolSplitter()

        db_data = '\n'.join(PEPTIDE_DB_CASE1)

        with io.BytesIO(db_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_database(handle)

        noncoding_data = ['>ENST0050|1', 'GAVPARRR']
        noncoding_fasta = '\n'.join(noncoding_data)

        with io.BytesIO(noncoding_fasta.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_database_noncoding(handle)

        self.assertIn('Noncoding', splitter.order.keys())
