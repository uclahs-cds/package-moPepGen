""" Module to test PeptideVariantGraph """
from typing import Tuple, Dict
import unittest
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import aa, dna, svgraph, vep


def create_pgraph(data:dict
        ) -> Tuple[svgraph.PeptideVariantGraph,Dict[int, svgraph.PeptideNode]]:
    """ Create a peptide variant graph from data """
    root = svgraph.PeptideNode(None)
    graph = svgraph.PeptideVariantGraph(root)
    node_list = {0: root}
    for key, val in data.items():
        seq = aa.AminoAcidSeqRecord(
            Seq(val[0]),
            id='ENST00001',
            transcript_id='ENST00001'
        )
        variants = []
        for it in val[2]:
            if it is None:
                continue
            location_transcirpt = FeatureLocation(start=it[0], end=it[1])
            location_peptide = FeatureLocation(start=it[5], end=it[6])
            var_record = vep.VariantRecord(
                    location_transcirpt,
                    ref=it[2],
                    alt=it[3]
                )
            var_record = vep.VEPVariantRecord(
                variant=var_record,
                consequences=it[4]
            )
            variant = svgraph.VariantRecordWithCoordinate(
                variant=var_record,
                location=location_peptide
            )
            variants.append(variant)
        node = svgraph.PeptideNode(seq, variants)
        node_list[key] = node
        for in_node_key in val[1]:
            node_list[in_node_key].add_out_edge(node)
    return graph, node_list


class TestPeptideVariantGraph(unittest.TestCase):    
    def test_expand_backward_no_inbound(self):
        """ > Test expand backward when the node has no inbound nodes.
        
                  V-P
                 /
            MNAAC-VC-PRN
                 \  /   
                  DC    
        """
        data = {
            1: ('MNAAC', [0], [None]),
            2: ('V', [1],[(0, 3, 'TCT', 'T', [], 0, 1)]),
            3: ('P', [2], [None]),
            4: ('VC', [1], [None]),
            5: ('DC', [1], [(0, 1, 'T', 'A', [], 0, 1)]),
            6: ('PRN', [4, 5], [None])
        }
        
        graph, nodes = create_pgraph(data)
        graph.rule = 'trypsin'
        graph.expand_alignment_backward(nodes[1])
        
        seqs = set([str(node.seq.seq) for node in graph.root.out_nodes])
        self.assertEqual(seqs, set(['MNAACV', 'MNAACVC', 'MNAACDC']))
        return
    
    def test_expand_backward_one_inbound(self):
        """ > Test expand backward when the node has no inbound nodes.
        
                       V-P
                      /
            MNAAC-GCVV-VC-PRN
                      \  /   
                       DC    
        """
        data = {
            1: ('MNAAC', [0], [None]),
            2: ('GCVV', [1], [None]),
            3: ('V', [2],[(0, 3, 'TCT', 'T', [], 0, 1)]),
            4: ('P', [3], [None]),
            5: ('VC', [2], [None]),
            6: ('DC', [2], [(0, 1, 'T', 'A', [], 0, 1)]),
            7: ('PRN', [5, 6], [None])
        }
        
        graph, nodes = create_pgraph(data)
        graph.rule = 'trypsin'
        graph.expand_alignment_backward(nodes[2])
        
        seqs = set([str(node.seq.seq) for node in nodes[1].out_nodes])
        self.assertEqual(seqs, set(['GCVVV', 'GCVVVC', 'GCVVDC']))
        return
    
    def test_expand_forward(self):
        """ > expand forward
                  GCVVV-P
                 /
            MNAAR-GCVVVC-PDNK-DCAGP
                 \      /   
                  GCVVDC    
        """
        data = {
            1: ('MNAAR', [0], [None]),
            2: ('GCVVV', [1], [(0, 3, 'TCT', 'T', [], 0, 1)]),
            3: ('P', [3], [None]),
            4: ('GCVVVC', [1], [None]),
            5: ('GCVVDC', [1], [(0, 1, 'T', 'A', [], 0, 1)]),
            6: ('PDNK', [4, 5], [None]),
            7: ('DCAGP', [6], [None])
        }
        graph, nodes = create_pgraph(data)
        graph.rule = 'trypsin'
        graph.expand_alignment_forward(nodes[6])
        seqs = set([str(node.seq.seq) for node in nodes[7].in_nodes])
        self.assertEqual(seqs, set(['GCVVVCPDNK', 'GCVVDCPDNK']))
        return
    
    def test_expand_forward_with_new_cleave_site(self):
        """ > expand forward when there is a new cleave site
                  GCVVV-P
                 /
            MNAAR-GCVVVC-PDNK-DCAGP
                 \      /   
                  GCVVRC    
        """
        data = {
            1: ('MNAAR', [0], [None]),
            2: ('GCVVV', [1], [(0, 3, 'TCT', 'T', [], 0, 1)]),
            3: ('P', [3], [None]),
            4: ('GCVVVC', [1], [None]),
            5: ('GCVVRC', [1], [(0, 1, 'T', 'A', [], 0, 1)]),
            6: ('PDNK', [4, 5], [None]),
            7: ('DCAGP', [6], [None])
        }
        graph, nodes = create_pgraph(data)
        graph.rule = 'trypsin'
        graph.expand_alignment_forward(nodes[6])
        seqs = set([str(node.seq.seq) for node in nodes[7].in_nodes])
        self.assertEqual(seqs, set(['GCVVVCPDNK', 'CPDNK']))
        seqs = set([str(node.seq.seq) for node in nodes[1].out_nodes])
        self.assertIn('GCVVR', seqs)
        return
    

if __name__ == '__main__':
    unittest.main()