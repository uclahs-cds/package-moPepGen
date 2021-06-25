""" Module for CircularVariantGraph """
import unittest
from Bio.Seq import Seq
from test import create_variant, create_dgraph2
from moPepGen import svgraph, seqvar, dna


class TestCVG(unittest.TestCase):
    """ Test Case for CVG """
    def test_create_cvg(self):
        """ Test the created cvg is cyclic """
        seq = Seq('AATTGGCCCCGGTTAA')
        locations = []
        seq = dna.DNASeqRecordWithCoordinates(seq, locations)
        graph = svgraph.CircularVariantGraph(seq, 'ENST0001')
        self.assertEqual(graph.root.seq.seq, seq.seq)

    def test_apply_variant(self):
        """ Apply variant to a CVG """
        seq = Seq('AATTGGCCCCGGTTAA')
        locations = []
        seq = dna.DNASeqRecordWithCoordinates(seq, locations)
        graph = svgraph.CircularVariantGraph(seq, 'ENST0001')
        data = (10, 11, 'G', 'T', 'SNV', '')
        variant = create_variant(*data)
        node = graph.root
        tail = graph.apply_variant(node, variant)
        self.assertEqual(str(tail.seq.seq), 'AATTGGCCCC')
        node = graph.root
        self.assertEqual(str(node.seq.seq), 'AATTGGCCCC')
        node = node.get_reference_next()
        self.assertEqual(str(node.seq.seq), 'G')
        node = node.get_reference_next()
        self.assertEqual(str(node.seq.seq), 'GTTAA')
        node = node.get_reference_next()
        self.assertIs(node, graph.root)

    def test_align_all_variants(self):
        r"""
                   A
                  / \
            ATAGGG-G-CCTGCT
        """
        data = {
            1: ('ATAGGG', [4], []),
            2: ('A', [1], [(0, 'G', 'A', 'SNV', '')]),
            3: ('G', [1], []),
            4: ('CCTGCT', [2, 3], [])
        }
        graph, _ = create_dgraph2(data, True)
        graph:svgraph.CircularVariantGraph
        graph.align_all_variants()
        self.assertEqual(str(graph.root.seq.seq), 'ATAGGG')

    def test_find_orf(self):
        r"""
                T
               / \
            ATA-G-GCTGCT
        """
        data = {
            1: ('ATA', [4], []),
            2: ('T', [1], [(0, 'G', 'T', 'SNV', '')]),
            3: ('G', [1], []),
            4: ('GCTGCT', [2, 3], [])
        }
        cvg, _ = create_dgraph2(data, True)
        cvg:svgraph.CircularVariantGraph
        cvg.attrs = {'id': 'ENSG0001-2'}
        tvg = cvg.find_all_orfs()
        self.assertIs(tvg.root.seq, None)
        self.assertEqual(len(tvg.root.out_edges), 1)
        for edge in tvg.root.out_edges:
            self.assertTrue(edge.out_node.seq.seq.startswith('A'))
