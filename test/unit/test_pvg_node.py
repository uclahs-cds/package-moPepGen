""" Module to test PVGNode """
from typing import Dict, List, Tuple
from collections import deque
import unittest
from moPepGen.svgraph import PVGNode
from test.unit import create_pvg_node, create_variant


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
        nodes = node.split_node_archipel(None)
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
        nodes = node.split_node_archipel(None)
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
        nodes = node.split_node_archipel(None)
        self.assertEqual(
            [str(x.seq.seq) for x in nodes],
            ['S', 'RS', 'K', 'Y', 'L', 'G']
        )

    def test_split_ref_to_single_amino_acid_4(self):
        """ Test split nodes with global variant """
        seq = 'SRSKYLG'
        variant_1 = (4, 5, 'A', '<INS>', 'Insertion', 'RI-100-200')
        variant_data = [
            (1, 7, *variant_1)
        ]
        node = create_pvg_node(seq, variants=variant_data)
        v = create_variant(*variant_1)
        nodes = node.split_node_archipel(v)
        self.assertEqual(
            [str(x.seq.seq) for x in nodes],
            ['S', 'R', 'S', 'K', 'Y', 'L', 'G']
        )

    def test_split_ref_to_single_amino_acid_5(self):
        """ Test split nodes with global variant and more """
        seq = 'SRSKYLG'
        variant_1 = (4, 5, 'A', '<INS>', 'Insertion', 'RI-100-200')
        variant_2 = (13, 14, 'A', 'AAAAAT', 'SNV', '')
        variant_data = [
            (1, 7, *variant_1),
            (1, 7, *variant_2)
        ]
        node = create_pvg_node(seq, variants=variant_data)
        v = create_variant(*variant_1)
        nodes = node.split_node_archipel(v)
        self.assertEqual(
            [str(x.seq.seq) for x in nodes],
            ['S', 'RSKYLG']
        )
