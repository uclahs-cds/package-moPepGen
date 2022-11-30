""" Test cases for TVGNodeCollapser """
import unittest
from test.unit import create_variant, create_dna_seq_with_coordinates
from moPepGen import seqvar
from moPepGen.svgraph import TVGNode
from moPepGen.svgraph.TVGNodeCollapser import TVGNodeCollapser
from moPepGen.SeqFeature import FeatureLocation


class TestTVGNodeCollapser(unittest.TestCase):
    """ Test cases for TVGNodeCollapser """
    def test_collapse_global_variant(self):
        """ This is to make sure that lobal variant is ignored. """
        var1 = create_variant(0, 100, 'A', '<circRNA>', 'circRNA', 'CIRC-0001')
        var2 = create_variant(4, 5, 'A', 'T', 'SNV', 'ENST0001-4-A-T')

        seq = create_dna_seq_with_coordinates('ACG')
        vars_coord_1 = [
            seqvar.VariantRecordWithCoordinate(var1, FeatureLocation(start=0, end=3)),
            seqvar.VariantRecordWithCoordinate(var2, FeatureLocation(start=1, end=2))
        ]
        node1 = TVGNode(seq, vars_coord_1, global_variant=var1)

        vars_coord_2 = [
            seqvar.VariantRecordWithCoordinate(var1, FeatureLocation(start=0, end=3))
        ]
        node2 = TVGNode(seq, vars_coord_2, global_variant=var1)

        collapser = TVGNodeCollapser()
        collapser.collapse(node1)
        res = collapser.collapse(node2)
        # node1 has 1 non global variant while node2 has 0, so node1 should be
        # discarded.
        self.assertIs(res, node1)
        self.assertIsNot(res, node2)

    def test_collapse_subgraphs(self):
        """ Tests that nodes in different subgraphs will not be collapsed. """
        var1 = create_variant(4, 5, 'A', 'T', 'SNV', 'ENST0001-4-A-T')
        seq = create_dna_seq_with_coordinates('ACG')
        vars_coord_1 = [
            seqvar.VariantRecordWithCoordinate(var1, FeatureLocation(start=0, end=1))
        ]

        node1 = TVGNode(seq, vars_coord_1, subgraph_id='abcde', level=1)
        node2 = TVGNode(seq, subgraph_id='ENST0001')

        collapser = TVGNodeCollapser()
        collapser.collapse(node1)
        res = collapser.collapse(node2)

        self.assertIsNone(res)
        self.assertEqual(len(collapser.pool), 2)
        self.assertIn(node1, set(collapser.mapper.values()))
        self.assertIn(node2, set(collapser.mapper.values()))
