""" Test module for the transcript graph """
import unittest
from typing import Deque
from collections import deque
from test import create_dgraph2, create_dgraph1
from moPepGen import svgraph, seqvar
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
        first_node = graph.root.get_reference_next()
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
        first_node = graph.root.get_reference_next()
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
        first_node = graph.root.get_reference_next()
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
        first_node = graph.root.get_reference_next()
        graph.align_variants(first_node)
        self.assertEqual(str(first_node.get_reference_next().seq.seq), 'TCT')
        variant_seqs = [str(edge.out_node.seq.seq) for edge \
            in first_node.out_edges]
        self.assertEqual(set(['T', 'TCT',]), set(variant_seqs))

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
        queue:Deque[svgraph.DNANode] = deque([first_node])
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
