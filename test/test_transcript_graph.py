""" Test module for the transcript graph """
import unittest
from typing import Deque
from collections import deque
from Bio.Seq import Seq
from test import create_dgraph2, create_dgraph1
from moPepGen import svgraph, seqvar, dna
from moPepGen.SeqFeature import FeatureLocation


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
            3: ['TGG', [1,2], []]
        }
        graph, nodes = create_dgraph2(data)
        graph.expand_alignments(nodes[0])
        variant_site_nodes = [edge.out_node for edge in nodes[0].out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertEqual(set(['ATG', 'AAT']), set(variant_site_seqs))

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
        cur, branches = graph.expand_alignments(nodes[1])
        variant_site_nodes = [edge.out_node for edge in nodes[1].out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertEqual(set(['GTCTGAC', 'GTG']), set(variant_site_seqs))
        self.assertIs(cur, nodes[5])
        self.assertEqual(len(branches), 1)
        self.assertEqual(str(branches[0].seq.seq), 'CCCT')

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
        cur, branches = graph.expand_alignments(nodes[1])
        self.assertEqual(str(nodes[4].seq.seq), 'CCT')
        self.assertEqual(len(branches), 0)
        self.assertIs(cur, nodes[4])

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
            2: ('AA', [0], [(0, 'AA', 'A', 'INDEL', '')], True),
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

    def test_find_overlaps(self):
        """ Correct farthest overlap node is found. """
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)
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
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (4, 5, 'T', 'A', 'SNV', '4:T-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root
        graph.align_variants(first_node)
        self.assertEqual(str(first_node.get_reference_next().seq.seq), 'TCT')
        variant_seqs = [str(edge.out_node.seq.seq) for edge \
            in first_node.out_edges]
        self.assertEqual(set(['T', 'TCT', 'ACT']), set(variant_seqs))

        # case 2
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', 'INDEL', '4:TCT-T'),
            (7, 8, 'G', 'A', 'SNV', '7:G-A')
        ]
        graph = create_dgraph1(seq, variants)
        first_node = graph.root
        graph.align_variants(first_node)
        self.assertEqual(str(first_node.get_reference_next().seq.seq), 'TCTG')
        variant_seqs = [str(edge.out_node.seq.seq) for edge \
            in first_node.out_edges]
        self.assertEqual(set(['T', 'TCTA', 'TCTG']), set(variant_seqs))

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

                 T--                     T-G-CCCT
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
        expected = {'T', 'TCTG', 'TCTA'}
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
        graph.fit_into_codons()
        pgraph = graph.translate()
        self.assertIsInstance(pgraph, svgraph.PeptideVariantGraph)


if __name__ == '__main__':
    unittest.main()
