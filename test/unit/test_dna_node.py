""" Test module for DNA Node """
import unittest
from test.unit import create_three_frame_tvg
from moPepGen import svgraph


class TestTVGNode(unittest.TestCase):
    """ Test case for DNA node """
    def test_deep_copy(self):
        """ Test that the deepcopy creates a copy of the node and its
        downstream nodes. """
        data = {
            1: ('ATGTGGC', ['RF0'], []),
            2: ('CA', [1], [(0, 'C', 'CA', 'INDEL', '')]),
            3: ('C', [1], []),
            4: ('TGCTG', [2,3], []),
            5: ('TCCGT', [4], [])
        }
        seq = 'ATGTGGCCTGCTGTCCGT'
        _, nodes = create_three_frame_tvg(data, seq)
        node_copy = nodes[2].deepcopy()
        self.assertEqual(node_copy.seq.seq, nodes[2].seq.seq)
        self.assertIsNot(node_copy, nodes[2])

        out_node:svgraph.TVGNode = list(node_copy.out_edges)[0].out_node
        self.assertEqual(out_node.seq.seq, nodes[4].seq.seq)
        self.assertIsNot(out_node, nodes[4])

        out_node:svgraph.TVGNode = list(out_node.out_edges)[0].out_node
        self.assertEqual(out_node.seq.seq, nodes[5].seq.seq)
        self.assertIsNot(out_node, nodes[5])

    def test_find_farthest_node_with_overlap_case1(self):
        r"""
                 T--
                /   \
            ATGG-TCT-G-CCCT
                    \ /
                     A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        seq = 'ATGGTCTGACCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[6])

    def test_find_farthest_node_with_overlap_case2(self):
        r""" Test case for the last node is too short. Should include the
        following node.
                 T--   T
                /   \ / \
            ATGG-TCT-G-C-C-CT
                    \ /
                     A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('T', [4], [(0, 'C', 'T', 'SNV', '')]),
            7: ('C', [4,5], []),
            8: ('C', [6,7], []),
            9: ('CT', [8], [])
        }
        seq = 'ATGGTCTGCCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        node = nodes[1].find_farthest_node_with_overlap()
        print(node.seq.seq)
        self.assertIs(node, nodes[9])

    def test_find_farthest_node_with_overlap_case3(self):
        r""" For a bubble that do not merge, should return None. This would  be
        the case of mutation at the last nucleotide.
                 C
                /
            ATGG-T
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], []),
            3: ['C', [1], [(0, 'T', 'C', 'SNV', '')]]
        }
        seq = 'ATGGT'
        _, nodes = create_three_frame_tvg(data, seq)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertEqual(str(node.seq.seq), 'T')


    def test_find_farthest_node_with_overlap_case5_with_branch(self):
        r"""
                 T-G-CCCT
                /
            ATGG-TCTGAC-G-CCCT
                       \ /
                        A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            5: ('TCTGAC', [1], []),
            6: ('G', [5], []),
            7: ('A', [5], [(0, 'G', 'T', 'SNV', '')]),
            8: ('CCCT', [6, 7], [])
        }
        seq = 'ATGGTCTGACGCCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[5])

    def test_find_farthest_node_with_exclusive_outbond(self):
        r"""
                 T--
                /    \
            ATGG-TCTC-G-CCCT-GTTGGCCC
                     \ /
                      A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            3: ('TCTC', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], []),
            7: ('GTTGGCCC', [5,6], []),
        }
        seq = 'ATGGTCTCGCCCTGTTGGCCC'
        _, nodes = create_three_frame_tvg(data, seq)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[7])
