""" Test module for DNA Node """
import unittest
from test.unit import create_three_frame_tvg, create_dna_seq_with_coordinates,\
    create_variant
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, svgraph


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

    def test_check_stop_altering_false(self):
        """ Test that non stop altering variants are not affected. """
        seq = create_dna_seq_with_coordinates(
            seq='AAAAAAACTTTT',
            locations=[((0,7),(0,7)), ((8,12), (8,12))],
            orf=[0, 12]
        )
        variant = create_variant(7,8,'G', 'C', 'SNV', '')
        variant = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=7, end=8)
        )
        node = svgraph.TVGNode(seq, [variant])
        node.check_stop_altering(12)
        self.assertFalse(variant.is_stop_altering)

    def test_check_stop_altering_true(self):
        """ Test that stop altering variants are identified. """
        seq = create_dna_seq_with_coordinates(
            seq='AAAAAAAC',
            locations=[((0,7),(0,7))],
            orf=[0, 7]
        )
        variant = create_variant(7,8,'G', 'C', 'SNV', '')
        variant = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=7, end=8)
        )
        node = svgraph.TVGNode(seq, [variant])
        node.check_stop_altering(7)
        self.assertTrue(variant.is_stop_altering)

    def test_check_stop_altering_true_noncoding(self):
        """ Test that stop altering variants are identified. """
        seq = create_dna_seq_with_coordinates(
            seq='AAAAAATAC',
            locations=[((0,8),(0,8))],
            orf=[0, 9]
        )
        variant = create_variant(8,9,'G', 'C', 'SNV', '')
        variant = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=8, end=9)
        )
        node = svgraph.TVGNode(seq, [variant])
        node.check_stop_altering(None)
        self.assertTrue(variant.is_stop_altering)
