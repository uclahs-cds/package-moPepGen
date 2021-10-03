""" Module for CircularVariantGraph """
from typing import List, Tuple
import unittest
from Bio.Seq import Seq
from moPepGen import svgraph, dna, circ
from moPepGen.SeqFeature import FeatureLocation


def create_circ_model(gene_id:str, fragments:List[Tuple[int,int]],
        _id:str, tx_ids:List[str]=None, gene_name:str=None, intron:List[int]=None
        ) -> circ.CircRNAModel:
    """ create circRNA model """
    fragments = [FeatureLocation(start=i, end=j) for i,j in fragments]
    return circ.CircRNAModel(
        gene_id=gene_id,
        fragments=fragments,
        intron=intron or [],
        _id=_id,
        transcript_ids=tx_ids or 'ENST0001',
        gene_name=gene_name or 'SYMB'
    )

class TestCVG(unittest.TestCase):
    """ Test Case for CVG """
    def test_init_three_frames(self):
        """ Test the created cvg is cyclic """
        circ_record = create_circ_model('ENSG0001', [(0,8),(10,18)], 'CIRCXXX')
        seq = Seq('AATTGGCCCCGGTTAA')
        locations = []
        seq = dna.DNASeqRecordWithCoordinates(seq, locations)
        graph = svgraph.ThreeFrameCVG(seq, 'ENST0001', circ_record=circ_record)
        graph.init_three_frames()
        for root in graph.reading_frames:
            self.assertTrue(root.get_reference_next().is_inbond_of(root))

    def test_extend_loop(self):
        """ extend loop """
        circ_record = create_circ_model('ENSG0001', [(0,8),(10,18)], 'CIRCXXX')
        seq = Seq('AATTGGCCCCGGTTAA')
        locations = []
        seq = dna.DNASeqRecordWithCoordinates(seq, locations)
        graph = svgraph.ThreeFrameCVG(seq, 'ENST0001', circ_record=circ_record)
        graph.init_three_frames()
        graph.extend_loop()
        for root in graph.reading_frames:
            node = list(list(root.out_edges)[0].out_node.out_edges)[0].out_node
            self.assertEqual(str(node.seq.seq), str(seq.seq))
