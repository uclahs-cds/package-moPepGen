""" Module for testing PeptidePoolSplitter """
import copy
import io
import unittest
from test.unit import create_aa_record
from moPepGen.err import VariantSourceNotFoundError
from moPepGen.aa import PeptidePoolSplitter, VariantSourceSet, \
    VariantPeptideInfo, VariantPeptidePool


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
    'ENST0001': {
        'SNV-1001-T-A': 'gSNP',
        'SNV-1002-T-A': 'gSNP',
        'SNV-1003-T-A': 'sSNV',
        'SNV-1004-T-A': 'sSNV',
        'INDEL-1101-TTTT-T': 'gINDEL',
        'INDEL-1102-TTTT-T': 'gINDEL',
        'INDEL-1103-TTTT-T': 'sINDEL',
        'INDEL-1104-TTTT-T': 'sINDEL'
    },
    'ENST0002': {
        'SNV-2001-T-A': 'gSNP',
        'SNV-2002-T-A': 'gSNP',
        'SNV-2003-T-A': 'sSNV',
        'SNV-2004-T-A': 'sSNV',
        'INDEL-2101-TTTT-T': 'gINDEL',
        'INDEL-2102-TTTT-T': 'gINDEL',
        'INDEL-2103-TTTT-T': 'sINDEL',
        'INDEL-2104-TTTT-T': 'sINDEL'
    },
    'ENST0003': {
        'SNV-3001-T-A': 'gSNP',
        'SNV-3002-T-A': 'gSNP',
        'SNV-3003-T-A': 'sSNV',
        'SNV-3004-T-A': 'sSNV',
        'INDEL-3101-TTTT-T': 'gINDEL',
        'INDEL-3102-TTTT-T': 'gINDEL',
        'INDEL-3103-TTTT-T': 'sINDEL',
        'INDEL-3104-TTTT-T': 'sINDEL'
    }
}

SOURCE_ORDER = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3}

class TestVariantSourceSet(unittest.TestCase):
    """ Test cases for VariantSourceSet """

    def setUp(self):
        """ Set up """
        super().setUp()
        VariantSourceSet.reset_levels()

    def test_source_levels(self):
        """ Test comparing variant source set """
        levels = copy.copy(SOURCE_ORDER)
        VariantSourceSet.set_levels(levels)
        set1 = VariantSourceSet(['gSNP'])
        self.assertEqual(set1.levels_map, levels)

        levels2 = copy.copy(levels)
        levels2['Fusion'] = 4
        set1.set_levels(levels2)
        set2 = VariantSourceSet(['sSNV'])
        self.assertEqual(set2.levels_map, levels2)

    def test_comparison(self):
        """ Test comparisons """
        levels = copy.copy(SOURCE_ORDER)
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
        levels = copy.copy(SOURCE_ORDER)
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

    def test_from_variant_peptide_case2(self):
        """ Test forming from variant peptide """
        levels = copy.copy(SOURCE_ORDER)
        VariantSourceSet.set_levels(levels)
        label_map = {'ENST0001': {'SNV-1157-G-A': 'sSNV'}}

        peptide = create_aa_record('KHIRJ','ENST0002|SNV-1157-G-A|1')
        with self.assertRaises(VariantSourceNotFoundError):
            VariantPeptideInfo.from_variant_peptide(peptide, label_map)

        peptide = create_aa_record('KHIRJ','ENST0001|INDEL-1120-GAAA-A|1')
        with self.assertRaises(VariantSourceNotFoundError):
            VariantPeptideInfo.from_variant_peptide(peptide, label_map)

class TestPeptidePoolSplitter(unittest.TestCase):
    """ Test cases for testing PeptidePoolSplitter """
    def test_append_order_noncoding(self):
        """ Test the group order noncoding """
        levels = copy.copy(SOURCE_ORDER)
        splitter = PeptidePoolSplitter(order=levels)
        splitter.append_order_noncoding()
        self.assertEqual(splitter.order['Noncoding'], 4)

    def test_load_gvf(self):
        """ test loading gvf """
        levels = copy.copy(SOURCE_ORDER)
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

    def test_split_database_case1(self):
        """ Test split database case 1. Single peptide, single source """
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=copy.copy(LABEL_MAP1),
            sources=copy.copy(list(SOURCE_ORDER.keys()))
        )
        splitter.split(1, [])
        self.assertIn('gSNP', splitter.databases.keys())
        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

    def test_split_database_case2(self):
        """ Test split database case 2 """
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=copy.copy(LABEL_MAP1),
            sources=copy.copy(list(SOURCE_ORDER.keys()))
        )
        splitter.split(1, [])

        remaining = PeptidePoolSplitter.get_remaining_database_key()
        self.assertEqual({'gSNP', remaining}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

        received = {str(x.seq) for x in splitter.databases[remaining].peptides}
        expected = {'SSSSSSSK'}
        self.assertEqual(expected, received)

    def test_split_database_case3(self):
        """ Test split database case 3. Two sources. """
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=copy.copy(LABEL_MAP1),
            sources=copy.copy(list(SOURCE_ORDER.keys()))
        )
        splitter.split(2, [])

        self.assertEqual({'gSNP', 'gSNP-sSNV'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

        received = {str(x.seq) for x in splitter.databases['gSNP-sSNV'].peptides}
        expected = {'SSSSSSSK'}
        self.assertEqual(expected, received)

    def test_split_database_case4(self):
        """ Test split database case 4. Two transcripts. """
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1 ENST0002|INDEL-2101-TTTT-T|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=copy.copy(LABEL_MAP1),
            sources=copy.copy(list(SOURCE_ORDER.keys()))
        )
        splitter.split(1, [])

        self.assertEqual({'gSNP', 'gINDEL'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

        peptide = list(splitter.databases['gINDEL'].peptides)[0]
        self.assertEqual('SSSSSSSK', str(peptide.seq))
        expected = 'ENST0002|INDEL-2101-TTTT-T|1 ENST0001|SNV-1002-T-A|SNV-1003-T-A|1'
        self.assertEqual(peptide.description, expected)

    def test_split_database_case5(self):
        """ Test split database case 5. Two transcripts with additional splitter. """
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [
                'SSSSSSSK',
                'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1'
                ' ENST0002|INDEL-2101-TTTT-T|INDEL-2104-TTTT-T|1'
            ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=copy.copy(LABEL_MAP1),
            sources=copy.copy(list(SOURCE_ORDER.keys()))
        )
        splitter.split(1, [{'sINDEL'}])

        self.assertEqual({'gSNP', 'sINDEL-additional'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

        peptide = list(splitter.databases['sINDEL-additional'].peptides)[0]
        self.assertEqual('SSSSSSSK', str(peptide.seq))
        expected = 'ENST0002|INDEL-2101-TTTT-T|INDEL-2104-TTTT-T|1' +\
            ' ENST0001|SNV-1002-T-A|SNV-1003-T-A|1'
        self.assertEqual(peptide.description, expected)
