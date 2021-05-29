from moPepGen.vep import VEPVariantRecord
import unittest
from collections import deque
from Bio.Seq import Seq
from moPepGen import svgraph
from moPepGen import vep
from moPepGen import dna
from moPepGen.SeqFeature import FeatureLocation


def create_graph(seq, variants) -> svgraph.TranscriptVariantGraph:
    location = dna.MatchedLocation(
        query=FeatureLocation(start=0, end=len(seq)),
        ref=FeatureLocation(start=0, end=len(seq))
    )
    seq = dna.DNASeqRecordWithCoordinates(
        seq=Seq(seq),
        locations=[location],
        orf=FeatureLocation(start=0, end=len(seq))
    )
    graph = svgraph.TranscriptVariantGraph(seq, '', None)
    records = []
    for start, end, ref, alt, consequences in variants:
        location = FeatureLocation(start=start, end=end)
        records.append(vep.VEPVariantRecord(
            variant=vep.VariantRecord(location=location, ref=ref, alt=alt),
            consequences=consequences
        ))
    graph.create_variant_graph(records)
    return graph
    

class TestTranscriptGraph(unittest.TestCase):
    
    def test_create_graph(self):
        seq = 'ATGGTCTGCCCTCTGAAC'
        variants = [
            (3, 4, 'G', 'T', ['missense']),
            (3, 4, 'G', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        variant_site_nodes = [edge.out_node for edge in graph.root.out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertIn('T', variant_site_seqs)
        self.assertIn('A', variant_site_seqs)
        self.assertIn('G', variant_site_seqs)
        
        # testing that the variant nodes all go to the same reference node.
        nodes = set([next(iter(edge.out_node.out_edges)).out_node for edge \
            in graph.root.out_edges])
        self.assertEqual(len(nodes), 1)

        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['missense']),
            (4, 5, 'T', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)

        return
        
    
    def test_expand_bubble(self):
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (3, 4, 'G', 'T', ['missense']),
            (3, 4, 'G', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        graph.expand_alignments(graph.root)
        variant_site_nodes = [edge.out_node for edge in graph.root.out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertEqual(set(['TTC', 'ATC', 'GTC']), set(variant_site_seqs))

        seqs = [edge.in_node.seq.seq for edge in graph.root\
            .get_reference_next().get_reference_next().in_edges]
        self.assertEqual(set(['TTC', 'ATC', 'GTC']), set(seqs))
        return
    
    def test_find_overlaps(self):
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['missense']),
            (4, 5, 'T', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        farthest = graph.root.find_farthest_node_with_overlap()
        self.assertEqual('GCCCTCTGAACTGA', str(farthest.seq.seq))

        # case 2
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['deletion']),
            (7, 8, 'G', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        farthest = graph.root.find_farthest_node_with_overlap()
        self.assertEqual('CCCTCTGAACTGA', str(farthest.seq.seq))
        return

    def test_align_variants(self):
        # case 1
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['missense']),
            (4, 5, 'T', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        graph.align_varints(graph.root)
        self.assertEqual(str(graph.root.get_reference_next().seq.seq), 'TCT')
        variant_seqs = [str(edge.out_node.seq.seq) for edge \
            in graph.root.out_edges]
        self.assertEqual(set(['T', 'TCT', 'ACT']), set(variant_seqs))

        # case 2
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['missense']),
            (7, 8, 'G', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        graph.align_varints(graph.root)
        self.assertEqual(str(graph.root.get_reference_next().seq.seq), 'TCTG')
        variant_seqs = [str(edge.out_node.seq.seq) for edge \
            in graph.root.out_edges]
        self.assertEqual(set(['TG', 'TCTG', 'TCTA']), set(variant_seqs))
        return
    
    def test_fit_into_codons(self):
        # case 1
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['missense']),
            (4, 5, 'T', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        graph.fit_into_codons()
        queue:deque[svgraph.Node] = deque([graph.root])
        while queue:
            cur = queue.pop()
            if cur.seq and cur.out_edges:
                self.assertEqual(len(cur.seq.seq) % 3, 0)
            for edge in cur.out_edges:
                queue.appendleft(edge.out_node)
        return

    def test_translate(self):
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['missense']),
            (4, 5, 'T', 'A', ['missense'])
        ]
        graph = create_graph(seq, variants)
        graph.fit_into_codons()
        pgraph = graph.translate()
        return

if __name__ == '__main__':
    unittest.main()