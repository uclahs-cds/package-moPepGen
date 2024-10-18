""" Module to test PVGNode """
from typing import Dict, List, Tuple
from collections import deque
import unittest
from moPepGen.svgraph import PVGNode
from test.unit import create_pvg_node


class TestPVGNode(unittest.TestCase):
    """ Test PVGNode """
    @staticmethod
    def collect_all_upstream_nodes(node:PVGNode):
        """ Collect all upstream nodes """
        cur = node
        nodes = deque([cur])
        while True:
            in_nodes = cur.get_in_nodes()
            if not in_nodes:
                return nodes
            cur = in_nodes[0]
            nodes.appendleft(cur)

    def test_split_ref_to_single_amino_acid_1(self):
        """ Test split to single amino acid """
        seq = 'SRSKYLG'
        node = create_pvg_node(seq)
        nodes = node.split_node_archipel()
        self.assertTrue(all(len(x.seq.seq) == 1 for x in nodes))
        self.assertEqual(''.join([str(x.seq.seq) for x in nodes]), seq)

    def test_split_ref_to_single_amino_acid_2(self):
        """ Test split to single amino acid with SNVs. All amino acids should be
        separated. """
        seq = 'SRSKYLG'
        variant_1 = (4, 5, 'A', 'T', 'SNV', '')
        variant_2 = (13, 14, 'A', 'T', 'SNV', '')
        variant_data = [
            (2, 3, *variant_1),
            (5, 6, *variant_2)
        ]
        node = create_pvg_node(seq, variants=variant_data)
        nodes = node.split_node_archipel()
        self.assertTrue(all(len(x.seq.seq) == 1 for x in nodes))
        self.assertEqual(''.join([str(x.seq.seq) for x in nodes]), seq)

    def test_split_ref_to_single_amino_acid_3(self):
        """ Test split to single amino acid with INDELs """
        seq = 'SRSKYLG'
        variant_1 = (4, 5, 'A', 'TT', 'INDEL', '')
        variant_2 = (13, 14, 'A', 'T', 'SNV', '')
        variant_data = [
            (1, 3, *variant_1),
            (5, 6, *variant_2)
        ]
        node = create_pvg_node(seq, variants=variant_data)
        nodes = node.split_node_archipel()
        self.assertEqual(
            [str(x.seq.seq) for x in nodes],
            ['S', 'RS', 'K', 'Y', 'L', 'G']
        )
