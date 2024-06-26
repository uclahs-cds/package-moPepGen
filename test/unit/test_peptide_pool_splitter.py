""" Module for testing PeptidePoolSplitter """
import copy
import io
import unittest
from test.unit import create_aa_record, create_genomic_annotation, get_tx2gene_and_coding_tx
from moPepGen.err import VariantSourceNotFoundError
from moPepGen.aa import PeptidePoolSplitter, VariantSourceSet, \
    VariantPeptideInfo, VariantPeptidePool
from moPepGen.aa.PeptidePoolSplitter import LabelSourceMapping


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
    'ENSG0001': {
        'SNV-1001-T-A': 'gSNP',
        'SNV-1002-T-A': 'gSNP',
        'SNV-1003-T-A': 'sSNV',
        'SNV-1004-T-A': 'sSNV',
        'INDEL-1101-TTTT-T': 'gINDEL',
        'INDEL-1102-TTTT-T': 'gINDEL',
        'INDEL-1103-TTTT-T': 'sINDEL',
        'INDEL-1104-TTTT-T': 'sINDEL',
        'FUSION-ENST0001:1050-ENST0003:3090': 'Fusion'
    },
    'ENSG0002': {
        'SNV-2001-T-A': 'gSNP',
        'SNV-2002-T-A': 'gSNP',
        'SNV-2003-T-A': 'sSNV',
        'SNV-2004-T-A': 'sSNV',
        'INDEL-2101-TTTT-T': 'gINDEL',
        'INDEL-2102-TTTT-T': 'gINDEL',
        'INDEL-2103-TTTT-T': 'sINDEL',
        'INDEL-2104-TTTT-T': 'sINDEL',
        'CIRC-ENST0002-E1-E2': 'circRNA',
        'SE-2100': 'altSplice'
    },
    'ENSG0003': {
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

SOURCE_ORDER = {'gSNP': 0, 'gINDEL': 1, 'sSNV': 2, 'sINDEL': 3, 'Fusion':4,
    'altSplice':5, 'circRNA': 6, 'NovelORF': 7}

ANNOTATION_ATTRS = [
    [
        {
            'gene_id': 'ENSG0001',
            'gene_type': 'protein_coding'
        }, {
            'transcript_id': 'ENST0001',
            'gene_id': 'ENSG0001',
            'protein_id': 'ENSP0001',
            'gene_type': 'protein_coding'
        }
    ], [
        {
            'gene_id': 'ENSG0002',
            'gene_type': 'protein_coding'
        }, {
            'transcript_id': 'ENST0002',
            'gene_id': 'ENSG0002',
            'protein_id': 'ENSP0002',
            'gene_type': 'protein_coding'
        }
    ], [
        {
            'gene_id': 'ENSG0003',
            'gene_type': 'protein_coding'
        }, {
            'transcript_id': 'ENST0003',
            'gene_id': 'ENSG0003',
            'protein_id': 'ENSP0003',
            'gene_type': 'protein_coding'
        }
    ], [
        {
            'gene_id': 'ENSG0004',
            'gene_type': 'protein_coding'
        }, {
            'transcript_id': 'ENST0004',
            'gene_id': 'ENSG0004',
            'protein_id': 'ENSP0004',
            'gene_type': 'protein_coding'
        }
    ]
]
ANNOTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0][0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 100, ANNOTATION_ATTRS[0][0]),
        'transcripts': ['ENST0001']
    },{
        'gene_id': ANNOTATION_ATTRS[1][0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (1000, 1100, ANNOTATION_ATTRS[0][0]),
        'transcripts': ['ENST0002']
    },{
        'gene_id': ANNOTATION_ATTRS[2][0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (2000, 2100, ANNOTATION_ATTRS[0][0]),
        'transcripts': ['ENST0003']
    },{
        'gene_id': ANNOTATION_ATTRS[3][0]['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (3000, 3100, ANNOTATION_ATTRS[0][0]),
        'transcripts': ['ENST0004']
    }],
    'transcripts': [{
        'transcript_id': ANNOTATION_ATTRS[0][1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (0, 100, ANNOTATION_ATTRS[0][1]),
        'exon': [(0, 100, ANNOTATION_ATTRS[0][1])],
        'cds': [(0, 100, ANNOTATION_ATTRS[0][1])]
    }, {
        'transcript_id': ANNOTATION_ATTRS[1][1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (1000, 1100, ANNOTATION_ATTRS[1][1]),
        'exon': [(1000, 1100, ANNOTATION_ATTRS[1][1])],
        'cds': [(1000, 1100, ANNOTATION_ATTRS[1][1])]
    }, {
        'transcript_id': ANNOTATION_ATTRS[2][1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (2000, 2100, ANNOTATION_ATTRS[2][1]),
        'exon': [(2000, 2100, ANNOTATION_ATTRS[2][1])],
        'cds': [(2000, 2100, ANNOTATION_ATTRS[2][1])]
    }, {
        'transcript_id': ANNOTATION_ATTRS[3][1]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (3000, 3100, ANNOTATION_ATTRS[3][1]),
        'exon': [(3000, 3100, ANNOTATION_ATTRS[3][1])],
        'cds': []
    }]
}

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

        levels[frozenset({'sINDEL', 'circRNA'})] = 7
        VariantSourceSet.set_levels(levels)
        set1 = VariantSourceSet(['gSNP', 'circRNA'])
        set2 = VariantSourceSet(['sINDEL', 'circRNA'])
        self.assertTrue(set1 > set2)

        levels = copy.copy(SOURCE_ORDER)
        levels.update({
            'NovelORF': 5,
            frozenset({'NovelORF', 'circRNA'}): 6,
            'circRNA': 7
        })
        VariantSourceSet.set_levels(levels)
        set1 = VariantSourceSet(['NovelORF', 'circRNA'])
        set2 = VariantSourceSet(['circRNA'])
        self.assertTrue(set1 < set2)

class TestVariantPeptideInfo(unittest.TestCase):
    """ Test VariantPeptideInfo """
    def test_from_variant_peptide(self):
        """ Test forming from variant peptide """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        levels = copy.copy(SOURCE_ORDER)
        VariantSourceSet.set_levels(levels)
        label_map_data = {'ENSG0001': {'SNV-1157-G-A': 'sSNV'}}
        label_map = LabelSourceMapping(label_map_data)
        peptide = create_aa_record('KHIRJ','ENST0001|SNV-1157-G-A|1')
        infos = VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)
        self.assertEqual([x.sources for x in infos], [{'sSNV'}])

        label_map_data = {
            'ENSG0001': {'SNV-1157-G-A': 'gSNP', 'INDEL-100-GGGGG-G': 'gINDEL'}
        }
        label_map = LabelSourceMapping(label_map_data)
        label = 'ENST0001|SNV-1157-G-A|INDEL-100-GGGGG-G|1'
        peptide = create_aa_record('KHIRJ', label)
        infos = VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)
        self.assertEqual([x.sources for x in infos], [{'gSNP', 'gINDEL'}])

        label_map_data = {
            'ENSG0001': {'SNV-1157-G-A': 'gSNP', 'INDEL-100-GGGGG-G': 'gINDEL'},
            'ENSG0002': {'SNV-254-T-A': 'sSNV'}
        }
        label_map = LabelSourceMapping(label_map_data)
        label = 'ENST0001|SNV-1157-G-A|INDEL-100-GGGGG-G|1 ENST0002|SNV-254-T-A|1'
        peptide = create_aa_record('KHIRJ', label)
        infos = VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)
        self.assertEqual([x.sources for x in infos], [{'gSNP', 'gINDEL'}, {'sSNV'}])

    def test_from_variant_peptide_case2(self):
        """ Test forming from variant peptide """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        levels = copy.copy(SOURCE_ORDER)
        VariantSourceSet.set_levels(levels)
        label_map_data = {'ENSG0001': {'SNV-1157-G-A': 'sSNV'}}
        label_map = LabelSourceMapping(label_map_data)

        peptide = create_aa_record('KHIRJ','ENST0002|SNV-1157-G-A|1')
        with self.assertRaises(VariantSourceNotFoundError):
            VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)

        peptide = create_aa_record('KHIRJ','ENST0001|INDEL-1120-GAAA-A|1')
        with self.assertRaises(VariantSourceNotFoundError):
            VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)

    def test_from_variant_peptide_noncoding(self):
        """ Test forming from noncidng peptide """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        levels = copy.copy(SOURCE_ORDER)
        VariantSourceSet.set_levels(levels)
        label_map_data = {'ENSG0001': {'SNV-1157-G-A': 'sSNV'}}
        label_map = LabelSourceMapping(label_map_data)

        peptide = create_aa_record('KHIRJ','ENST0004|ENSG0004|ORF1|1')
        infos = VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)
        self.assertIn('NovelORF', infos[0].sources)

        peptide = create_aa_record('KHIRJ','ENST0004|ORF1|1')
        infos = VariantPeptideInfo.from_variant_peptide(peptide, tx2gene, coding_tx, label_map)
        self.assertIn('NovelORF', infos[0].sources)

class TestPeptidePoolSplitter(unittest.TestCase):
    """ Test cases for testing PeptidePoolSplitter """
    def test_append_order_noncoding(self):
        """ Test the group order novel ORF """
        levels = copy.copy(SOURCE_ORDER)
        splitter = PeptidePoolSplitter(order=levels)
        splitter.append_order_internal_sources()
        self.assertEqual(splitter.order['NovelORF'], 7)

    def test_load_gvf(self):
        """ test loading gvf """
        levels = copy.copy(SOURCE_ORDER)
        splitter = PeptidePoolSplitter(order=levels)

        file_data = '\n'.join(['\t'.join(x) for x in GVF_CASE1])
        with io.BytesIO(file_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_gvf(handle)
        self.assertEqual(
            splitter.label_map.get_source('ENSG0001', 'SNV-611-G-A'), 'gSNP'
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
            splitter.label_map.get_source('ENSG0001', 'SNV-611-G-A'), 'gSNP'
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
        """ Test loading database of novel ORF """
        splitter = PeptidePoolSplitter()

        db_data = '\n'.join(PEPTIDE_DB_CASE1)

        with io.BytesIO(db_data.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_database(handle)

        noncoding_data = ['>ENST0050|1', 'GAVPARRR']
        noncoding_fasta = '\n'.join(noncoding_data)

        with io.BytesIO(noncoding_fasta.encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                splitter.load_database(handle)

        received = {str(x.seq) for x in splitter.peptides.peptides}
        self.assertIn('GAVPARRR', received)

    def test_split_database_case1(self):
        """ Test split database case 1. Single peptide, single source """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [[ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ]]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(1, [], tx2gene, coding_tx)
        self.assertIn('gSNP', splitter.databases.keys())
        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

    def test_split_database_case2(self):
        """ Test split database case 2 """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(1, [], tx2gene, coding_tx)

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
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(2, [], tx2gene, coding_tx)

        self.assertEqual({'gSNP', 'gSNP-sSNV'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

        received = {str(x.seq) for x in splitter.databases['gSNP-sSNV'].peptides}
        expected = {'SSSSSSSK'}
        self.assertEqual(expected, received)

    def test_split_database_case4(self):
        """ Test split database case 4. Two transcripts. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1 ENST0002|INDEL-2101-TTTT-T|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(1, [], tx2gene, coding_tx)

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
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            ['SSSSSSSK', 'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1'
             ' ENST0002|INDEL-2101-TTTT-T|INDEL-2104-TTTT-T|1']
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(1, [{'gSNP'}], tx2gene, coding_tx)

        self.assertEqual({'gSNP', 'gSNP-additional'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['gSNP'].peptides}
        expected = {'SSSSSSSR'}
        self.assertEqual(expected, received)

        peptide = list(splitter.databases['gSNP-additional'].peptides)[0]
        self.assertEqual('SSSSSSSK', str(peptide.seq))
        expected = (
            'ENST0001|SNV-1002-T-A|SNV-1003-T-A|1'
            + ' ENST0002|INDEL-2101-TTTT-T|INDEL-2104-TTTT-T|1'
        )
        self.assertEqual(peptide.description, expected)

    def test_split_database_circ_rna(self):
        """ Test split database with circRNA. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'CIRC-ENST0002-E1-E2|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(1, [{'circRNA'}], tx2gene, coding_tx)

        self.assertEqual({'gSNP', 'circRNA'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['circRNA'].peptides}
        expected = {'SSSSSSSK'}
        self.assertEqual(expected, received)

    def test_split_database_fusion(self):
        """ Test split database with fusion. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ],
            [ 'SSSSSSSK', 'FUSION-ENST0001:1050-ENST0003:3090|1' ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=copy.copy(SOURCE_ORDER),
            label_map=label_map
        )
        splitter.split(1, [{'Fusion'}], tx2gene, coding_tx)

        self.assertEqual({'gSNP', 'Fusion'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['Fusion'].peptides}
        expected = {'SSSSSSSK'}
        self.assertEqual(expected, received)

    def test_split_database_source_comb_order(self):
        """ Test split database with source order of combinations. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        anno.transcripts['ENST0005'] = copy.deepcopy(anno.transcripts['ENST0002'])
        anno.transcripts['ENST0005'].is_protein_coding = False
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [
                'SSSSSSSR',
                'CIRC-ENST0002-E1-E2|1 ENST0005|SE-2100|ORF2|1'
            ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        order = copy.copy(SOURCE_ORDER)
        order.update({
            'NovelORF': 6,
            frozenset(['altSplice', 'NovelORF']): 7,
            'circRNA': 8
        })
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=order,
            label_map=label_map
        )
        splitter.split(2, [], tx2gene, coding_tx)

        self.assertEqual({'altSplice-NovelORF'}, set(splitter.databases.keys()))

        received = {str(x.seq) for x in splitter.databases['altSplice-NovelORF'].peptides}
        expected = {x[0] for x in peptides_data}
        self.assertEqual(expected, received)

    def test_split_database_source_comb_order_case2(self):
        """ Test split database with source order of combinations of more than 2. """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        anno.transcripts['ENST0005'] = copy.deepcopy(anno.transcripts['ENST0002'])
        anno.transcripts['ENST0005'].is_protein_coding = False
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [
                'SSSSFSSR',
                'CIRC-ENST0002-E1-E2|1 ENST0005|SE-2100|W2F-5|ORF-2|1'
            ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        order = copy.copy(SOURCE_ORDER)
        order.update({
            'NovelORF': 6,
            frozenset(['altSplice', 'NovelORF']): 7,
            frozenset(['altSplice', 'NovelORF', 'CodonReassign']): 8,
            'circRNA': 9
        })
        splitter = PeptidePoolSplitter(
            peptides=peptides,
            order=order,
            label_map=label_map,
        )
        splitter.split(2, [{'altSplice', 'NovelORF'}], tx2gene, coding_tx)

        key = 'altSplice-NovelORF-additional'
        self.assertEqual({key}, set(splitter.databases.keys()))
        received = {str(x.seq) for x in splitter.databases[key].peptides}
        expected = {x[0] for x in peptides_data}
        self.assertEqual(expected, received)
