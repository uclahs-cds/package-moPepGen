""" Test module for the transcript graph """
import unittest
from typing import Deque
from collections import deque
from test.unit import create_dgraph2, create_dgraph1, create_dna_record_dict, \
    create_genomic_annotation, \
    create_variant, create_variants
from moPepGen import svgraph, seqvar
from moPepGen.SeqFeature import FeatureLocation

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

class TestTranscriptGraph(unittest.TestCase):
    """ Test case for the transcript graph """
    def test_splice(self):
        """ Test node splice at position """
        data = {
            1: ('AACCTTGG', [], [])
        }
        graph, nodes = create_dgraph2(data)
        left, right = graph.splice(nodes[1], 4, 'reference')
        self.assertIs(left, graph.root)
        self.assertEqual(str(left.seq.seq), 'AACC')
        self.assertIn(right, [e.out_node for e in left.out_edges])
        self.assertIn(left, [e.in_node for e in right.in_edges])

    def test_apply_variant_case1(self):
        """
                               C
                              / \
        AACCTTGG    ->    AACC-T-TGG
        """
        data = {
            1: ('AACCTTGG', [], [])
        }
        graph, nodes = create_dgraph2(data)
        variant = seqvar.VariantRecord(
            location=FeatureLocation(start=4, end=5),
            ref='T',
            alt='C',
            _type='SNV',
            _id=''
        )
        node = graph.apply_variant(nodes[1], variant)
        self.assertEqual(str(node.seq.seq), 'AACC')
        self.assertIs(node, graph.root)
        seqs = {str(e.out_node.seq.seq) for e in node.out_edges}
        seqs_expected = set(['T', 'C'])
        self.assertEqual(seqs, seqs_expected)

    def test_apply_fusion_case1(self):
        r""" Fusion breakpoint: acceptor 8, donor 8
             C                    A     A
            / \                  / \   / \
        AACC-T-TGGCGGTTC      TGG-G-TCC-T-TC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACC', [], []),
            2: ('T', [1], []),
            3: ('C', [1], [(0, 'C', 'T', 'SNV', '')]),
            4: ('TGGCGGTTC', [2,3], [])
        }
        graph, nodes = create_dgraph2(data)

        attrs = {
            'ACCEPTER_GENE_ID': 'ENSG0001',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST0001.1',
            'ACCEPTER_POSITION': 20,
            'ACCEPTER_CHROM': 'chr1'
        }
        var_fusion = create_variant(8, 9, 'C', '<FUSION>', 'Fusion',
            'FUSIONXXX', attrs=attrs)
        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (14, 15, 'G', 'A', 'SNV', '', None, 'ENST0001.1'),
            (18, 19, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        variants = {'ENST0001.1': create_variants(var_data)}
        variant_pool = seqvar.VariantRecordPool(transcriptional=variants)

        acceptor_tail = graph.apply_fusion(nodes[4], var_fusion,
            variant_pool, genome, anno)

        self.assertEqual(len(acceptor_tail.out_edges), 2)
        self.assertEqual(str(acceptor_tail.seq.seq), 'TGG')
        for edge in acceptor_tail.out_edges:
            if edge.type != 'reference':
                accepter_root = edge.out_node
        self.assertEqual(str(accepter_root.seq.seq), 'ATGG')
        self.assertEqual(len(accepter_root.in_edges), 1)

        # also test alignment will treat the first fusioned node as reference
        graph.align_variants(acceptor_tail)
        self.assertEqual(str(accepter_root.seq.seq), 'ATGG')

    def test_apply_insertion_case1(self):
        r"""
                                    CCCTATG
                                   /       \
        AACCTGGCGGTTC     ->    AAC-C-------TGGCGGTTC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACCTGGCGGTTC', [], [])
        }
        graph, nodes = create_dgraph2(data)

        attrs = {
            'GENE_ID': 'ENSG0001',
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        var_insertion = create_variant(3,4,'C','<INS>','Insertion','', attrs)

        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (19, 20, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        variants = {'ENST0001.1': create_variants(var_data)}
        variant_pool = seqvar.VariantRecordPool(transcriptional=variants)

        node = graph.apply_insertion(nodes[1], var_insertion,
            variant_pool, genome, anno)
        self.assertEqual(str(node.seq.seq), 'AAC')
        node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        self.assertEqual(node_seqs, {'C', 'CCCTATG'})


    def test_apply_insertion_case2(self):
        r"""
                                       T    A
                                      / \  / \
                                    CC-C-TA-T-G
                                   /           \
        AACCTGGCGGTTC     ->    AAC-C-----------TGGCGGTTC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACCTGGCGGTTC', [], [])
        }
        graph, nodes = create_dgraph2(data)

        attrs = {
            'GENE_ID': 'ENSG0001',
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        var_insertion = create_variant(3,4,'C','<INS>','Insertion','', attrs)

        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (8, 9, 'C', 'T', 'SNV', '', None, 'ENST0001.1'),
            (11, 12, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        variants = {'ENST0001.1': create_variants(var_data)}
        variant_pool = seqvar.VariantRecordPool(transcriptional=variants)

        node = graph.apply_insertion(nodes[1], var_insertion,
            variant_pool, genome, anno)

        self.assertEqual(str(node.seq.seq), 'AAC')
        node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        self.assertEqual(node_seqs, {'C', 'CC'})

        node = list(filter(lambda x:str(x.out_node.seq.seq) == 'CC', node.out_edges))[0].out_node
        node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        self.assertEqual(node_seqs, {'C', 'T'})

        node = list(filter(lambda x:str(x.out_node.seq.seq) == 'C', node.out_edges))[0].out_node
        node = list(node.out_edges)[0].out_node
        node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        self.assertEqual(node_seqs, {'T', 'A'})

    def test_apply_substitution_case1(self):
        r"""                         T
                                    / \
                                   C-C-TATG
                                  /         \
        AACCTGGCGGTTC    ->    AAC-CTGGC-----GGTTC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACCTGGCGGTTC', [], [])
        }
        graph, nodes = create_dgraph2(data)

        attrs = {
            'GENE_ID': 'ENSG0001',
            'START': 3,
            'END': 8,
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        var_sub = create_variant(3,8,'C','<SUB>','Substitution','', attrs)

        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (8, 9, 'C', 'T', 'SNV', '', None, 'ENST0001.1'),
            (19, 20, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }

        variants = {'ENST0001.1': create_variants(var_data)}
        variant_pool = seqvar.VariantRecordPool(transcriptional=variants)

        node = graph.apply_substitution(nodes[1], var_sub, variant_pool,
            genome, anno)

        self.assertEqual(str(node.seq.seq), 'AAC')
        node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        self.assertEqual(node_seqs, {'CTGGC', 'C'})

        node = list(filter(lambda x:str(x.out_node.seq.seq) == 'C', node.out_edges))[0].out_node
        node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        self.assertEqual(node_seqs, {'C', 'T'})


    def test_apply_deletion_case1(self):
        r"""
                                   C-------
                                  /        \
        AACCTGGCGGTTC    ->    AAC-CTGGC----GGTTC
        """
        seq = 'AACCTGGCGGTTC'
        var_data = [
            (3,8,'C','<DEL>','Deletion','')
        ]
        varaints = create_variants(var_data)
        del_var = varaints.pop(-1)
        graph = create_dgraph1(seq, [])
        graph.add_null_root()
        node = graph.root.get_reference_next()
        node = graph.apply_variant(node, del_var)
        self.assertEqual(str(node.seq.seq), 'AAC')
        self.assertEqual(len(node.out_edges), 2)
        for edge in node.out_edges:
            x = edge.out_node
            if edge.type == 'reference':
                self.assertEqual(str(x.seq.seq), 'CTGGC')
            else:
                self.assertEqual(str(x.seq.seq), 'C')
                self.assertEqual(len(x.in_edges), 1)
                self.assertEqual(len(x.out_edges), 1)

    def test_create_dgraph1(self):
        """ Transcript grpah is created """
        seq = 'ATGGTCTGCCCTCTGAAC'
        variants = [
            (3, 4, 'G', 'T', 'SNV', '3:G-T'),
            (3, 4, 'G', 'A', 'SNV', '3:G-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root
        variant_site_nodes = [edge.out_node for edge in first_node.out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertIn('T', variant_site_seqs)
        self.assertIn('A', variant_site_seqs)
        self.assertIn('G', variant_site_seqs)

        # testing that the variant nodes all go to the same reference node.
        nodes = {next(iter(edge.out_node.out_edges)).out_node for edge \
            in first_node.out_edges}
        self.assertEqual(len(nodes), 1)

        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)

    def test_truncate_left_variant_right(self):
        """ test truncate left when variant on the right"""
        data = {
            1: ('ATGG', [], []),
            2: ('CCATG', [1], [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1], []),
            4: ('CCCT', [2,3], [])
        }
        _, nodes = create_dgraph2(data)
        left = nodes[2].truncate_left(2)
        self.assertEqual(str(left.seq.seq),'CC')
        self.assertEqual(len(left.variants), 0)
        self.assertEqual(str(nodes[2].seq.seq),'ATG')
        self.assertEqual(len(nodes[2].variants), 1)
        self.assertEqual(nodes[2].variants[0].location.start, 0)
        self.assertTrue(nodes[2].is_inbond_of(nodes[4]))
        self.assertFalse(left.in_edges)
        self.assertFalse(left.out_edges)

    def test_truncate_left_variant_left(self):
        """ test truncate left when variant on the left """
        data = {
            1: ('ATGG', [], []),
            2: ('CCATG', [1], [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1], []),
            4: ('CCCT', [2,3], [])
        }
        _, nodes = create_dgraph2(data)
        left = nodes[2].truncate_left(3)
        self.assertEqual(str(left.seq.seq),'CCA')
        self.assertEqual(len(left.variants), 1)
        self.assertEqual(left.variants[0].location.start, 2)
        self.assertEqual(str(nodes[2].seq.seq),'TG')
        self.assertEqual(len(nodes[2].variants), 0)

    def test_truncate_right_variant_left(self):
        """ test truncate left when variant on the left """
        data = {
            1: ('ATGG', [], []),
            2: ('CCATG', [1], [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1], []),
            4: ('CCCT', [2,3], [])
        }
        _, nodes = create_dgraph2(data)
        right = nodes[2].truncate_right(3)
        self.assertEqual(str(nodes[2].seq.seq),'CCA')
        self.assertEqual(len(nodes[2].variants), 1)
        self.assertEqual(nodes[2].variants[0].location.start, 2)
        self.assertEqual(str(right.seq.seq),'TG')
        self.assertEqual(len(right.variants), 0)
        self.assertTrue(nodes[1].is_inbond_of(nodes[2]))
        self.assertFalse(right.in_edges)
        self.assertFalse(right.out_edges)

    def test_truncate_right_variant_right(self):
        """ test truncate left when variant on the right """
        data = {
            1: ('ATGG', [], []),
            2: ('CCATG', [1], [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1], []),
            4: ('CCCT', [2,3], [])
        }
        _, nodes = create_dgraph2(data)
        right = nodes[2].truncate_right(2)
        self.assertEqual(str(nodes[2].seq.seq),'CC')
        self.assertEqual(len(nodes[2].variants), 0)
        self.assertEqual(str(right.seq.seq),'ATG')
        self.assertEqual(len(right.variants), 1)
        self.assertEqual(right.variants[0].location.start, 0)

    def test_expand_alignments_case1(self):
        """ Variant alignments are expanded """
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (3, 4, 'G', 'T', 'SNV', '3:G-T'),
            (3, 4, 'G', 'A', 'SNV', '3:G-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root
        graph.expand_alignments(first_node)
        variant_site_nodes = [edge.out_node for edge in first_node.out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertEqual(set(['TTC', 'ATC', 'GTC']), set(variant_site_seqs))

        seqs = [edge.in_node.seq.seq for edge in first_node\
            .get_reference_next().get_reference_next().in_edges]
        self.assertEqual(set(['TTC', 'ATC', 'GTC']), set(seqs))

    def test_expand_alignments_case2(self):
        """ When there is a vairant at the first nucleotide
                 AA-TGG
                /
            Null-A-TGG
        """
        data = {
            0: (None, [], []),
            1: ('A', [0], []),
            2: ('AA', [0], [(0, 'AA', 'A', 'INDEL', '')], True),
            3: ['TGG', [1], []],
            4: ['TGG', [2], []]
        }
        graph, nodes = create_dgraph2(data)
        graph.expand_alignments(nodes[0])
        variant_site_nodes = [edge.out_node for edge in nodes[0].out_edges]
        variant_site_seqs = {str(node.seq.seq) for node in variant_site_nodes}
        self.assertEqual({'ATGG', 'AATGG'}, variant_site_seqs)

    def test_merge_with_outbounds(self):
        """
            ATGG-T-G-CCCT
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], []),
            3: ('G', [2], []),
            4: ('CCCT', [3], [])
        }
        graph, nodes = create_dgraph2(data)
        returned_nodes = graph.merge_with_outbonds(nodes[3])
        self.assertEqual(str(nodes[4].seq.seq), 'GCCCT')
        self.assertTrue(nodes[2].is_inbond_of(nodes[4]))
        self.assertFalse(nodes[2].is_inbond_of(nodes[3]))
        self.assertFalse(nodes[3].is_inbond_of(nodes[4]))
        self.assertEqual(returned_nodes[0], nodes[4])

    def test_expand_alignments_case3(self):
        r"""
                 T-G-CCCT
                /
            ATGG-TCTGAC-G-CCCT
                       \ /
                        A
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            3: ('G', [2], []),
            4: ('CCCT', [3], []),
            5: ('TCTGAC', [1], []),
            6: ('G', [5], []),
            7: ('A', [5], [(0, 'G', 'T', 'SNV', '')]),
            8: ('CCCT', [6, 7], [])
        }
        graph, nodes = create_dgraph2(data)
        end_nodes = graph.expand_alignments(nodes[1])
        variant_site_nodes = [edge.out_node for edge in nodes[1].out_edges]
        variant_site_seqs = {str(node.seq.seq) for node in variant_site_nodes}
        self.assertEqual({'GTCTGAC', 'GTGCCCT'}, variant_site_seqs)
        self.assertEqual(len(end_nodes), 2)
        seqs = {str(node.seq.seq) for node in end_nodes}
        expected = {'GTCTGAC', 'GTGCCCT'}
        self.assertEqual(seqs, expected)

    def test_expand_alignments_case4(self):
        r""" Tests that the node after variants are truncated.
                 A
                / \
            ATGG-G-CCCT
        """
        data = {
            1: ('ATGG', [], []),
            2: ('A', [1], [(0, 'G', 'A', 'SNV', '')]),
            3: ('G', [1], []),
            4: ('CCCT', [2, 3], [])
        }
        graph, nodes = create_dgraph2(data)
        end_nodes = graph.expand_alignments(nodes[1])
        self.assertEqual(str(nodes[4].seq.seq), 'CCT')
        self.assertEqual(len(end_nodes), 1)
        self.assertIn(nodes[4], end_nodes)

    def test_expand_alignments_case6(self):
        r"""
                 TCATGTA           In this test case, there are three groups
                /       \          of nodes outbond to node 1 (ATGG). Node
               | TCATGTG-CCCT      2 and 3 are siblings, so the bubble should
               |/                  expand and return the trailing CCCT (4)
               |  TG-CCCT          Node 7 and 8 are also siblings do they
               |/                  should expand and return CCT (9). Node 5
            ATGG-TCGTG-CCCT        is alone, so it should merge 6 and return
                \     /            the merged node (TGCCCT)
                 TCGTA
        """
        data = {
            1: ('ATGG', [], []),
            2: ('TCGTG', [1], []),
            3: ('TCGTA', [1], [(4, 'G', 'A', 'SNV', '')]),
            4: ('CCCT', [2,3], []),
            5: ('TG', [1], [(0, 'TCGT', 'T', 'INDEL', '')], True),
            6: ('CCCT', [5], []),
            7: ('TCATGTG', [1], [(1, 'C', 'CAT', 'INDEL', '')], True),
            8: ('TCATGTA', [1], [(1, 'C', 'CAT', 'INDEL', ''), (6, 'G', 'A', 'SNV', '')], True),
            9: ('CCCT', [7,8], [])
        }
        graph, nodes = create_dgraph2(data)
        end_nodes = graph.expand_alignments(nodes[1])

        variant_site_nodes = [edge.out_node for edge in nodes[1].out_edges]
        variant_site_seqs = {str(node.seq.seq) for node in variant_site_nodes}
        expected = {'GTCGTA', 'GTCGTG', 'GTGCCCT', 'GTCATGTGC', 'GTCATGTAC'}
        self.assertEqual(variant_site_seqs, expected)

        end_node_seqs = {str(node.seq.seq) for node in end_nodes}
        expected = {'CCCT', 'CCT', 'GTGCCCT'}
        self.assertEqual(end_node_seqs, expected)

    def test_find_orf_known_case1(self):
        r""" Test when start codon at 0 and no start lost mutation
            ATGG-G-CCCT
                \ /
                 A
        """
        data = {
            0: ('ATGG', [], []),
            1: ('A', [0], [(0, 'G', 'A', 'SNV', '')]),
            2: ('G', [0], []),
            3: ['CCCT', [1,2], []]
        }
        graph, nodes = create_dgraph2(data)
        graph.add_null_root()
        graph.seq.orf = FeatureLocation(start=0, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[0]))

    def test_find_orf_known_case2(self):
        r""" When start codon is not at at 0 and no start lost variants

            GGATGG-G-CCCT
                  \ /
                   A
        """
        data = {
            0: ('GGATGG', [], []),
            1: ('A', [0], [(0, 'G', 'A', 'SNV', '')]),
            2: ('G', [0], []),
            3: ['CCCT', [1,2], []]
        }
        graph, nodes = create_dgraph2(data)
        graph.add_null_root()
        graph.seq.orf = FeatureLocation(start=2, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[0]))
        self.assertTrue(nodes[0].seq.seq.startswith('ATG'))

    def test_find_orf_known_case3(self):
        r""" Variants before start codon
              T
             / \
            G-G-ATGG-G-CCCT
                    \ /
                     A
        """
        data = {
            0: ('G', [], []),
            1: ('T', [0], [(0, 'T', 'T', 'SNV', '')]),
            2: ('G', [0], []),
            3: ('ATGG', [1,2], []),
            4: ('A', [3], [(0, 'G', 'A', 'SNV', '')]),
            5: ('G', [3], []),
            6: ('CCCT', [4,5], [])
        }
        graph, nodes = create_dgraph2(data)
        graph.add_null_root()
        graph.seq.orf = FeatureLocation(start=2, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[3]))

    def test_find_orf_known_case4(self):
        r""" Start lost variants
                G
               / \
            GGA-T-GG-G-CCCT
                    \ /
                     A
        """
        data = {
            1: ('GGA', [], []),
            2: ('G', [1], [(0, 'T', 'G', 'SNV', '')]),
            3: ('T', [1], []),
            4: ('GG', [2,3], []),
            5: ('A', [4], [(0, 'G', 'A', 'SNV', '')]),
            6: ('G', [4], []),
            7: ('CCCT', [5,6], [])
        }
        graph, nodes = create_dgraph2(data)
        graph.add_null_root()
        graph.seq.orf = FeatureLocation(start=2, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[4]))
        self.assertTrue(nodes[4].seq.seq.startswith('ATG'))

    def test_find_orf_known_case5(self):
        r""" Start lost variants
                 AA
                /  \
            Null-A--TGG
        """
        data = {
            0: (None, [], []),
            1: ('A', [0], []),
            2: ('AA', [0], [(0, 'AA', 'A', 'INDEL', '')], True),
            3: ['TGG', [1,2], []]
        }
        graph, nodes = create_dgraph2(data)
        graph.seq.orf = FeatureLocation(start=0, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[3]))
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))

    def test_find_orf_known_case6(self):
        r""" Deletion at beginning
                 CG--
                /    \
            Null-ATGG-TGG
        """
        data = {
            0: (None, [], []),
            1: ('ATGG', [0], []),
            2: ('CG', [0], [(0, 'ATGG', 'CG', 'INDEL', '')], True),
            3: ['TGG', [1,2], []]
        }
        graph, nodes = create_dgraph2(data)
        graph.seq.orf = FeatureLocation(start=0, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[1]))
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))

    def test_find_orf_known_case7(self):
        r""" Deletion at beginning but cds start unknown
                 T---
                /    \
            Null-ATGG-TGG
        """
        data = {
            0: (None, [], []),
            1: ('ATGG', [0], []),
            2: ('T', [0], [(0, 'ATGG', 'T', 'INDEL', '')], True),
            3: ['TGG', [1,2], []]
        }
        graph, nodes = create_dgraph2(data, cds_start_nf=True)
        graph.seq.orf = FeatureLocation(start=0, end=len(graph.seq))
        graph.find_orf_known()
        self.assertTrue(graph.root.is_inbond_of(nodes[1]))
        seqs = {str(edge.out_node.seq.seq) for edge in graph.root.out_edges}
        expected = {'TTGG', 'ATGG'}
        self.assertEqual(seqs, expected)

    def test_find_orf_unknown_case1(self):
        r""" Test when start codon at 0 and no start lost mutation
            ATGG-G-CCCT
                \ /
                 A
        """
        data = {
            0: ('ATGG', [], []),
            1: ('A', [0], [(0, 'G', 'A', 'SNV', '')]),
            2: ('G', [0], []),
            3: ['CCCT', [1,2], []]
        }
        graph, _ = create_dgraph2(data)
        graph.add_null_root()
        graph.find_orf_unknown()
        self.assertEqual(len(graph.root.out_edges), 1)
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))

    def test_find_orf_unknown_case2(self):
        r""" When start codon is not at at 0 and no start lost variants

            GGATGG-G-CCCT
                  \ /
                   A
        """
        data = {
            0: ('GGATGG', [], []),
            1: ('A', [0], [(0, 'G', 'A', 'SNV', '')]),
            2: ('G', [0], []),
            3: ['CCCT', [1,2], []]
        }
        graph, _ = create_dgraph2(data)
        graph.add_null_root()
        graph.find_orf_unknown()
        self.assertEqual(len(graph.root.out_edges), 1)
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))

    def test_find_orf_unknown_case3(self):
        r""" Variants before start codon
              T
             / \
            G-G-ATGG-G-CCCT
                    \ /
                     A
        """
        data = {
            0: ('G', [], []),
            1: ('T', [0], [(0, 'T', 'T', 'SNV', '')]),
            2: ('G', [0], []),
            3: ('ATGG', [1,2], []),
            4: ('A', [3], [(0, 'G', 'A', 'SNV', '')]),
            5: ('G', [3], []),
            6: ('CCCT', [4,5], [])
        }
        graph, _ = create_dgraph2(data)
        graph.add_null_root()
        graph.find_orf_unknown()
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))

    def test_find_orf_unknown_case4(self):
        r"""
                 AA
                /  \
            Null-A--TGG
        """
        data = {
            0: (None, [], []),
            1: ('A', [0], []),
            2: ('AA', [0], [(0, 'AA', 'A', 'INDEL', '')]),
            3: ['TGG', [1,2], []]
        }
        graph, _ = create_dgraph2(data)
        graph.seq.orf = FeatureLocation(start=0, end=len(graph.seq))
        graph.find_orf_unknown()
        self.assertEqual(len(graph.root.out_edges), 2)
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))

    def test_find_orf_unknown_case5(self):
        r""" Variants before start codon

            GGATGG-G-TGCT
                  \ /
                   A
        """
        data = {
            0: ('GGATGG', [], []),
            1: ('A', [0], [(0, 'G', 'A', 'SNV', '')]),
            2: ('G', [0], []),
            3: ('TGCT', [1,2], [])
        }
        graph, _ = create_dgraph2(data)
        graph.add_null_root()
        graph.find_orf_unknown()
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))
        self.assertEqual(len(graph.root.out_edges), 2)

    def test_find_orf_unkonwn_frameshift(self):
        r""" finding unkonwn ORF with framefhift mutations

            GGCTGG-CGGTATGC-TGCT
                  \         /
                   C--------
        """
        data = {
            0: ('GGCTGG', [], []),
            1: ('C', [0], [(0, 'CGGTATGC', 'C', 'INDEL', '')]),
            2: ('CGGTATGC', [0], []),
            3: ('TGCT', [1,2], [])
        }
        graph, _ = create_dgraph2(data)
        graph.add_null_root()
        graph.find_orf_unknown()
        for edge in graph.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('ATG'))
        self.assertEqual(len(graph.root.out_edges), 1)

    def test_find_orf_unkonwn_after_mutation(self):
        """ ORF after mutation """
        var_1 = (0, 'CGG', 'C', 'INDEL', '', 0, 1)
        var_2 = (0, 'ACC', 'A', 'INDEL', '', 6, 7)
        data = {
            0: ('GGCTGG', [], []),
            1: ('CTCTGCA', [0], [var_1, var_2]),
            2: ('CGGTCTGCACC', [0], []),
            3: ('TGCT', [1,2], [])
        }
        graph,_ = create_dgraph2(data)
        graph.add_null_root()
        graph.find_orf_unknown()
        self.assertEqual(len(graph.root.out_edges), 1)
        orf = list(graph.root.out_edges)[0].out_node
        self.assertEqual(len(orf.variants), 1)
        self.assertEqual(list(orf.variants)[0].location.start, 0)

    def test_find_overlaps(self):
        r""" Correct farthest overlap node is found.
             T--
            /    \
        ATGG-T-CT-GCCCTCTGAACTGA
            \ /
             A
        """
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)
        graph.add_null_root()
        first_node = graph.root.get_reference_next()
        farthest = first_node.find_farthest_node_with_overlap()
        self.assertEqual('GCCCTCTGAACTGA', str(farthest.seq.seq))

        # case 2
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (7, 8, 'G', 'A', 'SNV', '7:G-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root.get_reference_next()
        farthest = first_node.find_farthest_node_with_overlap()
        self.assertEqual('CCCTCTGAACTGA', str(farthest.seq.seq))

    def test_align_variants_case1(self):
        """ Variants are algined. """
        # case 1
        # case 2
        #      T--
        #     /   \
        # ATGG-T-CT-GCCCTCTGAACTGA
        #     \ /
        #      A
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root
        graph.align_variants(first_node)
        self.assertEqual(str(first_node.get_reference_next().seq.seq), 'TCT')
        variant_seqs = {str(edge.out_node.seq.seq) for edge \
            in first_node.out_edges}
        self.assertEqual({'T', 'TCT', 'ACT'}, variant_seqs)

        # case 2
        #      T--
        #     /   \
        # ATGG-TCT-G-CCCTCTGAACTGA
        #         \ /
        #          A
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (7, 8, 'G', 'A', 'SNV', '7:G-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root
        graph.align_variants(first_node)
        self.assertEqual(str(first_node.get_reference_next().seq.seq), 'TCTG')
        variant_seqs = {str(edge.out_node.seq.seq) for edge \
            in first_node.out_edges}
        self.assertEqual({'TG', 'TCTA', 'TCTG'}, variant_seqs)

    def test_align_variants_case2(self):
        """ When there is a mutation at the first nucleotide.
                 AA
                /  \
            Null-A--TGG
        """
        data = {
            0: (None, [], []),
            1: ('A', [0], []),
            2: ('AA', [0], [(0, 'AA', 'A', 'INDEL', '')]),
            3: ['TGG', [1,2], []]
        }
        graph, nodes = create_dgraph2(data)
        node = graph.align_variants(nodes[0])
        self.assertIs(graph.root, node)
        self.assertIs(nodes[0], node)

    def test_align_variants_case3(self):
        r""" Start lost variants

                 T--                     TG-CCCT
                /   \                   /
            ATGG-TCT-G-CCCT   ->    ATGG-TCTG-CCCT
                    \ /                 \    /
                     A                   TCTA
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('A', [3], [(0, 'G', 'A', 'SNV', '')]),
            5: ('G', [2,3], []),
            6: ('CCCT', [4,5], [])
        }
        graph, nodes = create_dgraph2(data)
        node = graph.align_variants(nodes[1])
        seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        expected = {'TG', 'TCTG', 'TCTA'}
        self.assertEqual(seqs, expected)

    def test_align_variants_case4(self):
        r"""
                                               TCATGA
                                              /      \
                                             | TCATGT-GCCCT
                                             |/
                                             | TGT
                                             |/   \
                 T----                       | TGA-GCCCT
                /     \                      |/
            ATGG-T-C---G-T-GCCCT    ->    ATGG-TCGT-GCCCT
                  \   / \ /                   \    /
                   CAT   A                     TCGA
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCG', 'T', 'INDEL', '')]),
            3: ('T', [1], []),
            4: ('CAT', [3], [(0, 'C', 'CAT', 'INDEL', '')]),
            5: ('C', [3], []),
            6: ('G', [2,4,5], []),
            7: ('T', [6], []),
            8: ('A', [6], [(0, 'A', 'A', 'SNV', '')]),
            9: ('GCCCT', [7,8], [])
        }
        graph, nodes = create_dgraph2(data)
        node = graph.align_variants(nodes[1])
        seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
        expected = {'TCGA', 'TCGT', 'TGA', 'TGT', 'TCATGT', 'TCATGA'}
        self.assertEqual(seqs, expected)

    def test_fit_into_codons(self):
        """ Transcript variatn graph is fit into codons. """
        # case 1
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)
        graph.fit_into_codons()
        first_node = graph.root.get_reference_next()
        queue:Deque[svgraph.TVGNode] = deque([first_node])
        while queue:
            cur = queue.pop()
            if cur.seq and cur.out_edges:
                self.assertEqual(len(cur.seq.seq) % 3, 0)
            for edge in cur.out_edges:
                queue.appendleft(edge.out_node)

    def test_translate(self):
        """ Test the transcript graph is translated to a peptide graph """
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)
        graph.add_null_root()
        graph.find_orf_known()
        graph.fit_into_codons()
        pgraph = graph.translate()
        self.assertIsInstance(pgraph, svgraph.PeptideVariantGraph)
