""""""
import unittest
from test import create_dgraph


class TestDNANode(unittest.TestCase):
    def test_create_graph(self):
        data = {
            1: ('ATGTGGC', [], []),
            2: ('A', [1], [(0, 'C', 'A', 'SNV', '')]),
            3: ('C', [1], []),
            4: ('TGCTG', [2,3], [])
        }
        graph, nodes = create_dgraph(data)
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
            4: ('TGCTG', [2,3], [])
        }
        _, nodes = create_dgraph(data)
        node_copy = nodes[2].deepcopy()
        self.assertEqual(node_copy.seq.seq, nodes[2].seq.seq)
        self.assertIsNot(node_copy, nodes[2])
        out_node = next(iter(node_copy.out_edges)).out_node
        self.assertEqual(out_node.seq.seq, nodes[4].seq.seq)
        self.assertIsNot(out_node, nodes[4])
    
    def test_find_farthest_node_with_overlap_case1(self):
        """
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
        _, nodes = create_dgraph(data)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, nodes[6])

    def test_find_farthest_node_with_overlap_case2(self):
        """ Test case for the last node is too short. Should include the
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
        _, nodes = create_dgraph(data)
        node = nodes[1].find_farthest_node_with_overlap()
        print(node.seq.seq)
        self.assertIs(node, nodes[9])
    
    def test_find_farthest_node_with_overlap_case3(self):
        """ For a bubble that do not merge, should return None. This would  be
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
        graph, nodes = create_dgraph(data)
        node = nodes[1].find_farthest_node_with_overlap()
        self.assertIs(node, None)

if __name__ == '__main__':
    unittest.main()