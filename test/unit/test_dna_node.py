""" Test module for DNA Node """
import unittest
from test.unit import create_dgraph2
from moPepGen import svgraph


class TestTVGNode(unittest.TestCase):
    """ Test case for DNA node """
    def test_create_graph(self):
        """ Test the graph can be constructed successfully """
        data = {
            1: ('ATGTGGC', [], []),
            2: ('A', [1], [(0, 'C', 'A', 'SNV', '')]),
            3: ('C', [1], []),
            4: ('TGCTG', [2,3], [])
        }
        graph, nodes = create_dgraph2(data)
        self.assertEqual(str(graph.seq.seq), 'ATGTGGCCTGCTG')
        self.assertEqual(nodes[1].seq.seq, 'ATGTGGC')
        node = next(iter(nodes[3].out_edges)).out_node
        self.assertEqual(node.seq.seq, nodes[4].seq.seq)

    def test_deep_copy(self):
        """ Test that the deepcopy creates a copy of the node and its
        downstream nodes. """
        data = {
            1: ('ATGTGGC', [], []),
            2: ('CA', [1], [(0, 'C', 'CA', 'INDEL', '')]),
            3: ('C', [1], []),
            4: ('TGCTG', [2,3], []),
            5: ('TCCGT', [4], [])
        }
        _, nodes = create_dgraph2(data)
        node_copy = nodes[2].deepcopy()
        self.assertEqual(node_copy.seq.seq, nodes[2].seq.seq)
        self.assertIsNot(node_copy, nodes[2])

        out_node:svgraph.TVGNode = list(node_copy.out_edges)[0].out_node
        self.assertEqual(out_node.seq.seq, nodes[4].seq.seq)
        self.assertIsNot(out_node, nodes[4])
        self.assertEqual(len(out_node.frameshifts), 1)

        out_node:svgraph.TVGNode = list(out_node.out_edges)[0].out_node
        self.assertEqual(out_node.seq.seq, nodes[5].seq.seq)
        self.assertIsNot(out_node, nodes[5])
        self.assertEqual(len(out_node.frameshifts), 1)

    def test_next_node_to_branch_out_cas1(self):
        r""" Returns the frameshifting mutation
                 T--
                /   \
            ATGG-TCT-G-CCCT
                    \ /
                     A
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        _, nodes = create_dgraph2(data)
        node = nodes[1].next_node_to_branch_out(nodes[6])
        self.assertIs(node, nodes[2])

    def test_next_node_to_branch_out_case2(self):
        r""" Should return None for in frame mutation
                 T--
                /   \
            ATGG-TCGT-G-CCCT
                     \ /
                      A
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCGT', 'T', 'INDEL', '')]),
            3: ('TCGT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        _, nodes = create_dgraph2(data)
        node = nodes[1].next_node_to_branch_out(nodes[6])
        self.assertIs(node, None)

    def test_next_node_to_branch_out_case3(self):
        r""" Should return for long mutation even if it's inframe
                 T--
                /   \
            ATGG-TCGGGGGGGT-G-CCCT
                           \ /
                            A
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCGGGGGGGT', 'T', 'INDEL', '')]),
            3: ('TCGGGGGGGT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        _, nodes = create_dgraph2(data)
        node = nodes[1].next_node_to_branch_out(nodes[6], 5)
        self.assertIs(node, nodes[2])

    def test_find_farthest_node_with_overlap_case1(self):
        r"""
                 T--
                /   \
            ATGG-TCT-G-CCCT
                    \ /
                     A
        """
        data = {
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        _, nodes = create_dgraph2(data)
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
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('T', [4], [(0, 'C', 'T', 'SNV', '')]),
            7: ('C', [4,5], []),
            8: ('C', [6,7], []),
            9: ('CT', [8], [])
        }
        _, nodes = create_dgraph2(data)
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
            1: ('ATGG', [], []),
            2: ('T', [1], []),
            3: ['C', [1], [(0, 'T', 'C', 'SNV', '')]]
        }
        _, nodes = create_dgraph2(data)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertEqual(str(node.seq.seq), 'T')

    def test_find_farthest_node_with_overlap_case4_null_root(self):
        r""" For mutation at the first nucleotide.
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
        _, nodes = create_dgraph2(data)
        node = nodes[0].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[3])

    def test_find_farthest_node_with_overlap_case5_with_branch(self):
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
        _, nodes = create_dgraph2(data)
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
            1: ('ATGG', [], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            3: ('TCTC', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], []),
            7: ('GTTGGCCC', [5,6], []),
        }
        _, nodes = create_dgraph2(data)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[7])

    def test_find_farthest_node_with_circular(self):
        r"""
                 T--
                /    \
            $-ATGG-TCTC-G-CCCT-^
                     \ /
                      A
        """
        data = {
            1: ('ATGG', [6], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            3: ('TCTC', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        _, nodes = create_dgraph2(data, True)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[6])
