""" Module to test PeptideVariantGraph """
from typing import Tuple, Dict
import unittest
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import aa, svgraph, seqvar


def create_pgraph(data:dict, _id:str
        ) -> Tuple[svgraph.PeptideVariantGraph,Dict[int, svgraph.PVGNode]]:
    """ Create a peptide variant graph from data """
    root = svgraph.PVGNode(None)
    graph = svgraph.PeptideVariantGraph(root, _id)
    node_list = {0: root}
    for key, val in data.items():
        seq = aa.AminoAcidSeqRecord(
            Seq(val[0]),
            _id='ENST00001',
            transcript_id='ENST00001'
        )
        variants = []
        frameshifts = set()
        for it in val[2]:
            if it is None:
                continue
            location_transcirpt = FeatureLocation(start=it[0], end=it[1])
            location_peptide = FeatureLocation(start=it[6], end=it[7])
            var_record = seqvar.VariantRecord(
                location=location_transcirpt,
                ref=it[2],
                alt=it[3],
                _type=it[4],
                _id=it[5]
            )
            variant = seqvar.VariantRecordWithCoordinate(
                variant=var_record,
                location=location_peptide
            )
            variants.append(variant)
            if it[8]:
                frameshifts.add(variant.variant)
        node = svgraph.PVGNode(seq, variants, frameshifts=frameshifts)
        if len(val) >= 4:
            node.orf = val[3]
        node_list[key] = node
        for in_node_key in val[1]:
            node_list[in_node_key].add_out_edge(node)
    return graph, node_list


class TestPeptideVariantGraph(unittest.TestCase):
    """ Test case for peptide variant graph """
    def test_expand_backward_no_inbound(self):
        r""" > Test expand backward when the node has no inbound nodes.

                  V-P
                 /
            MNAAC-VC-PRN
                 \  /
                  DC
        """
        data = {
            1: ('MNAAC', [0], [None]),
            2: ('V', [1],[(0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)]),
            3: ('P', [2], [None]),
            4: ('VC', [1], [None]),
            5: ('DC', [1], [(0, 1, 'T', 'A', 'SNV', '0:T-A', 0, 1, False)]),
            6: ('PRN', [4, 5], [None])
        }

        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        graph.expand_alignment_backward(nodes[1])

        seqs = {str(node.seq.seq) for node in graph.root.out_nodes}
        self.assertEqual(seqs, {'MNAACV', 'MNAACVC', 'MNAACDC'})
        self.assertEqual(len(nodes[1].in_nodes), 0)
        self.assertEqual(len(nodes[1].out_nodes), 0)
        self.assertNotIn(nodes[1], graph.root.out_nodes)

    def test_expand_backward_one_inbound(self):
        r""" > Test expand backward when the node has no inbound nodes.

                       V-P
                      /
            MNAAC-GCVV-VC-PRN
                      \  /
                       DC
        """
        data = {
            1: ('MNAAC', [0], [None]),
            2: ('GCVV', [1], [None]),
            3: ('V', [2],[(0,3,'TCT','T','INDEL', '0:TCT-T', 0, 1, True)]),
            4: ('P', [3], [None]),
            5: ('VC', [2], [None]),
            6: ('DC', [2], [(0, 1, 'T', 'A', 'SNV', '0:T-A', 0, 1, False)]),
            7: ('PRN', [5, 6], [None])
        }

        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches = graph.expand_alignment_backward(nodes[2])

        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        self.assertEqual(seqs, {'GCVVV', 'GCVVVC', 'GCVVDC'})
        self.assertEqual(len(nodes[2].in_nodes), 0)
        self.assertEqual(len(nodes[2].out_nodes), 0)
        self.assertNotIn(nodes[2], nodes[1].out_nodes)

        self.assertEqual(len(branches), 2)
        seqs = {str(branch.seq.seq) for branch in branches}
        self.assertEqual(seqs, {'GCVVV', 'PRN'})

    def test_expand_forward(self):
        r""" > expand forward
                  GCVVV-P
                 /
            MNAAR-GCVVVC-PDNK-DCAGP
                 \      /
                  GCVVDC
        """
        data = {
            1: ('MNAAR', [0], [None]),
            2: ('GCVVV', [1], [(0, 3, 'TCT', 'T', 'INDEL', '', 0, 1, True)]),
            3: ('P', [3], [None]),
            4: ('GCVVVC', [1], [None]),
            5: ('GCVVDC', [1], [(0, 1, 'T', 'A', 'SNV', '', 0, 1, False)]),
            6: ('PDNK', [4, 5], [None]),
            7: ('DCAGP', [6], [None])
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        graph.expand_alignment_forward(nodes[6])
        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        self.assertEqual(seqs, {'GCVVVCPDNK', 'GCVVDCPDNK'})
        self.assertEqual(len(nodes[6].in_nodes), 0)
        self.assertEqual(len(nodes[6].out_nodes), 0)
        self.assertNotIn(nodes[6], nodes[4].out_nodes)
        self.assertNotIn(nodes[6], nodes[5].out_nodes)
        self.assertNotIn(nodes[6], nodes[7].in_nodes)

    def test_expand_forward_with_new_cleave_site(self):
        r""" > expand forward when there is a new cleave site
                  GCVVV-P
                 /
            MNAAR-GCVVVC-PDNK-DCAGP
                 \      /
                  GCVVRC
        """
        data = {
            1: ('MNAAR', [0], [None]),
            2: ('GCVVV', [1], [(0, 3, 'TCT', 'T', 'INDEL', '', 0, 1, True)]),
            3: ('P', [3], [None]),
            4: ('GCVVVC', [1], [None]),
            5: ('GCVVRC', [1], [(0, 1, 'T', 'A', 'SNV', '', 0, 1, False)]),
            6: ('PDNK', [4, 5], [None]),
            7: ('DCAGP', [6], [None])
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        graph.expand_alignment_forward(nodes[6])
        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        self.assertEqual(seqs, {'GCVVVCPDNK', 'CPDNK'})
        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        self.assertIn('GCVVR', seqs)
        self.assertEqual(len(nodes[6].in_nodes), 0)
        self.assertEqual(len(nodes[6].out_nodes), 0)
        self.assertNotIn(nodes[6], nodes[4].out_nodes)
        self.assertNotIn(nodes[6], nodes[5].out_nodes)
        self.assertNotIn(nodes[6], nodes[7].in_nodes)

    def test_merge_join(self):
        r"""
                                       NCWHSTQV
                                      /        \
            NCWV     V               | NCWVSTQV |
           /    \   / \              |/        \|
        AER-NCWH-STQ-Q-PK    ->    AER-NCWHSTQQ-PK
                                      \        /
                                       NCWVSTQQ
        """
        data = {
            1: ('AER', [0], [None]),
            2: ('NCWH', [1], [None]),
            3: ('NCWV', [1], [(0, 1, 'A', 'T', 'SNV', '', 0, 1, False)]),
            4: ('STQ', [2,3], []),
            5: ('Q', [4], []),
            6: ('V', [4], [(0, 1, 'A', 'T', 'SNV', '', 0, 1, False)]),
            7: ('PK', [5,6], [])
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches = graph.merge_join_alignments(nodes[4])
        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        seqs_expect = {'NCWHSTQV', 'NCWVSTQV', 'NCWHSTQQ', 'NCWVSTQQ'}
        self.assertEqual(seqs, seqs_expect)
        self.assertEqual(len(nodes[4].in_nodes), 0)
        self.assertEqual(len(nodes[4].out_nodes), 0)
        seqs = {str(node.seq.seq) for node in branches}
        expected = {'PK'}
        self.assertEqual(seqs, expected)

    def test_merge_join_with_cleave(self):
        r"""
                                       NCWHSTQV-
                                      /         \
            NCWR     V               | NCWR-STQV |
           /    \   / \              |/         \|
        AER-NCWH-STQ-Q-PK    ->    AER-NCWHSTQQ--PK
                                      \         /
                                       NCWR-STQQ
        """
        data = {
            1: ('AER', [0], [None]),
            2: ('NCWH', [1], [None]),
            3: ('NCWR', [1], [(0, 1, 'A', 'T', 'SNV', '', 0, 1, False)]),
            4: ('STQ', [2,3], []),
            5: ('Q', [4], []),
            6: ('V', [4], [(0, 1, 'A', 'T', 'SNV', '', 0, 1, False)]),
            7: ('PK', [5,6], [])
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches = graph.merge_join_alignments(nodes[4])
        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        seqs_expect = {'NCWHSTQV', 'NCWR', 'NCWHSTQQ', 'NCWR'}
        self.assertEqual(seqs, seqs_expect)
        self.assertEqual(len(nodes[4].in_nodes), 0)
        self.assertEqual(len(nodes[4].out_nodes), 0)
        seqs = {str(node.seq.seq) for node in branches}
        expected = {'PK'}
        self.assertEqual(seqs, expected)
        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        seqs_expect = {'NCWHSTQV', 'STQV', 'NCWHSTQQ', 'STQQ'}
        self.assertEqual(seqs, seqs_expect)

    def test_cross_join(self):
        r"""
            NCWV           V                 NCWVSTEEK-LPAQV
           /    \         / \               /         X     \
        AER-NCWH-STEEKLPAQ-Q-PK    ->    AER-NCWHSTEEK-LPAQQ-PK
        """
        data = {
            1: ('AER', [0], [None]),
            2: ('NCWH', [1], [None]),
            3: ('NCWV', [1], [(0, 1, 'A', 'T', 'SNV', '', 3, 4, False)]),
            4: ('STEEKLPAQ', [2,3], [None]),
            5: ('Q', [4], [None]),
            6: ('V', [4], [(0, 1, 'A', 'T', 'SNV', '', 0, 1, False)]),
            7: ('PK', [5,6], [None])
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches = graph.cross_join_alignments(nodes[4], 5)

        seqs = {str(node.seq.seq) for node in branches}
        expected = {'PK'}
        self.assertEqual(seqs, expected)

        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        seqs_expected = {'NCWVSTEEK', 'NCWHSTEEK'}
        self.assertEqual(seqs, seqs_expected)

        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        seqs_expected = {'LPAQV', 'LPAQQ'}
        self.assertEqual(seqs, seqs_expected)

        seqs = {str(node.seq.seq) for node in nodes[2].out_nodes}
        seqs_expected = {'LPAQV', 'LPAQQ'}
        self.assertEqual(seqs, seqs_expected)

        seqs = {str(node.seq.seq) for node in nodes[3].out_nodes}
        seqs_expected = {'LPAQV', 'LPAQQ'}
        self.assertEqual(seqs, seqs_expected)

    def test_call_variant_peptides_ambiguous_amino_acid(self):
        """ test call peptides with ambiguous amino acid """
        data = {
            1: ('SSSSK', [0], [None]),
            2: ('SSSSR', [1],[(0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)]),
            3: ('SSSIR', [1], [None]),
            4: ('SSXPK', [2,3], [None])
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_miscleavages(self):
        """ test micleavages is handled correctly """
        data = {
            1: ('SSSSK', [0], [None]),
            2: ('SSSSR', [1],[(0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)]),
            3: ('SSSIR', [1], [None]),
            4: ('SSSPK', [2,3], [None]),
            5: ('SSSVK', [4], [None])
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR', 'SSSSRSSSPK'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_truncated(self):
        """ test tuncated peptide is handled properly """
        data = {
            1: ('SSSSK', [0], [None]),
            2: ('SSSSR', [1],[(0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)]),
            3: ('SSSIR', [1], [None]),
            4: ('SSSPK', [2,3], [None]),
            5: ('SSSVK', [4], [None])
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        nodes[5].truncated = True
        peptides = graph.call_variant_peptides(2)
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR', 'SSSSKSSSSRSSSPK', 'SSSSRSSSPK'}
        self.assertEqual(received, expected)

    def test_call_peptides_check_orf(self):
        """ test calling peptides when checking for ORF """
        data = {
            1: ('SSSSK', [0], [None], [1,None]),
            2: ('SSSSJ', [0], [None], [2,None]),
            3: ('SSSSR', [1], [None], [1,None]),
            4: ('SSSSR', [2], [None], [2,None])
        }
        graph_id = 'ENST0001'
        graph, _ = create_pgraph(data, graph_id)
        peptides = graph.call_variant_peptides(0, check_variants=False, check_orf=True)

        received = {str(x.seq) for x in peptides}
        expected = {'SSSSK', 'SSSSR', 'SSSSJ'}
        self.assertEqual(received, expected)

        expected = [
            {f'{graph_id}|ORF1|1'}, {f'{graph_id}|ORF2|1'},
            {f'{graph_id}|ORF1|2', f'{graph_id}|ORF2|2'}
        ]
        received = [set(str(x.description).split(' ')) for x in peptides]
        self.assertEqual(len(expected), len(received))
        for it in received:
            self.assertIn(it, expected)
        for it in expected:
            self.assertIn(it, received)

    def test_call_peptides_id(self):
        """ test peptide ids are named correctly. """
        data = {
            1: ('SSSSK', [0], [None]),
            2: ('SSSSR', [1],[(0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)]),
            3: ('SSSIR', [1], [None]),
            4: ('SSXPK', [2,3], [None])
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        for peptide in peptides:
            self.assertTrue(peptide.description.startswith(graph.id))
