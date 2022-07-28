""" Module for CircularVariantGraph """
from typing import List, Tuple
import unittest
from Bio.Seq import Seq
from moPepGen import svgraph, dna, circ
from moPepGen.SeqFeature import FeatureLocation


def create_circ_model(tx_id:str, fragments:List[Tuple[int,int]],
        _id:str, gene_id:str=None, gene_name:str=None, intron:List[int]=None
        ) -> circ.CircRNAModel:
    """ create circRNA model """
    fragments = [FeatureLocation(start=i, end=j) for i,j in fragments]
    return circ.CircRNAModel(
        transcript_id=tx_id,
        fragments=fragments,
        intron=intron or [],
        _id=_id,
        gene_id=gene_id or 'ENSG0001',
        gene_name=gene_name or 'SYMB'
    )

class TestCVG(unittest.TestCase):
    """ Test Case for CVG """
    def test_init_three_frames(self):
        """ Test the created cvg is cyclic """
        circ_record = create_circ_model('ENST0001', [(0,8),(10,18)], 'CIRCXXX')
        seq = Seq('AATTGGCCCCGGTTAA')
        locations = []
        seq = dna.DNASeqRecordWithCoordinates(seq, locations)
        graph = svgraph.ThreeFrameCVG(seq, 'ENST0001', circ_record=circ_record)
        graph.init_three_frames()
        for root in graph.reading_frames:
            self.assertTrue(root.get_reference_next().is_inbond_of(root))

    def test_extend_loop(self):
        """ Extend loop. This test case ensures that the ThreeFrameCVG is
        extended into the linear structure by copying the graph for three
        times. """
        circ_record = create_circ_model('ENST0001', [(0,8),(10,18)], 'CIRCXXX')
        seq = Seq('AATTGGCCCCGGTTAA')
        locations = []
        seq = dna.DNASeqRecordWithCoordinates(seq, locations)
        graph = svgraph.ThreeFrameCVG(seq, 'ENST0001', circ_record=circ_record)
        graph.init_three_frames()
        graph.extend_loop()
        for cur in graph.reading_frames:
            n_loops = 0
            while cur.out_edges:
                self.assertEqual(len(cur.out_edges), 1)
                cur = cur.get_out_nodes()[0]
                self.assertEqual(cur.seq.seq, seq.seq)
                n_loops += 1
            self.assertEqual(n_loops, 4)
