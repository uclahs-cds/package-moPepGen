""" Module to test PeptideVariantGraph """
from typing import Tuple, Dict, List
import unittest
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen import aa, svgraph, seqvar


def create_pgraph(data:dict, _id:str, has_known_orf:bool=True
        ) -> Tuple[svgraph.PeptideVariantGraph,Dict[int, svgraph.PVGNode]]:
    """ Create a peptide variant graph from data """
    root = svgraph.PVGNode(None)
    graph = svgraph.PeptideVariantGraph(root, _id, has_known_orf)
    node_list = {0: root}
    for key, val in data.items():
        locs = []
        for (query_start, query_end), (ref_start, ref_end) in val[3]:
            loc = MatchedLocation(
                query=FeatureLocation(start=query_start, end=query_end),
                ref=FeatureLocation(start=ref_start, end=ref_end)
            )
            locs.append(loc)

        if len(val) >= 5:
            orf = val[4]
        else:
            orf = [0, None]

        seq = aa.AminoAcidSeqRecordWithCoordinates(
            Seq(val[0]),
            _id='ENST00001',
            transcript_id='ENST00001',
            locations=locs,
            orf=orf
        )
        seq.__class__ = aa.AminoAcidSeqRecordWithCoordinates
        variants:List[seqvar.VariantRecordWithCoordinate] = []
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
        node.orf = orf
        node_list[key] = node
        for in_node_key in val[1]:
            node_list[in_node_key].add_out_edge(node)

    return graph, node_list


class TestPeptideVariantGraph(unittest.TestCase):
    """ Test case for peptide variant graph """
    def test_find_reference_next(self):
        """ Test that the reference next and prev can be found """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        data = {
            1: ('SSSR', [0], [None], [((0,4), (0,4))]),
            2: ('SSSV', [1], [variant_1], [((0,3), (4,7))]),
            3: ('SSSG', [1], [None], [((0,4), (4,8))]),
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[1].find_reference_next(), nodes[3])
        graph.add_stop(nodes[2])
        self.assertIs(nodes[2].find_reference_next(), graph.stop)

    def test_find_reference_next_chimeric(self):
        """ When nodes are chimeric, such as combined because of frameshifting
        mutations. """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (5, 8, 'ACC', 'A', 'INDEL', '0:ACC-A', 0, 1, True)
        data = {
            1: ('SSSR', [0], [variant_1], [((1,4),(1,4))]),
            2: ('SSSV', [1],[variant_1, variant_2], [((1,4),(1,4))]),
            3: ('SSSG', [1], [variant_1], [((1,4),(1,4))]),
        }
        _, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[1].find_reference_next(), nodes[3])

    def test_find_reference_prev(self):
        """ Test that the reference next and prev can be found """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        data = {
            1: ('SSSR', [0], [None], [((0,4),(0,4))]),
            2: ('SSSV', [0],[variant_1], [((0,4),(0,4))]),
            3: ('SSSG', [1,2], [None], [((0,4),(0,4))]),
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[3].find_reference_prev(), nodes[1])
        self.assertIs(nodes[2].find_reference_prev(), graph.root)

    def test_find_reference_prev_chimeric(self):
        """ When nodes are chimeric, such as combined because of frameshifting
        mutations. """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (5, 8, 'ACC', 'A', 'INDEL', '0:ACC-A', 0, 1, True)
        data = {
            1: ('SSSR', [0], [variant_1], [((0,4),(0,4))]),
            2: ('SSSV', [0],[variant_1, variant_2], [((0,4),(0,4))]),
            3: ('SSSG', [1,2], [variant_1], [((0,4),(0,4))]),
        }
        _, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[3].find_reference_prev(), nodes[1])

    def test_expand_backward_no_inbound(self):
        r""" > Test expand backward when the node has no inbound nodes.

                  V-P
                 /
            MNAAC-VC-PRN
                 \  /
                  DC
        """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAC', [0], [None], locations),
            2: ('V', [1],[variant_1], locations),
            3: ('P', [2], [None], locations),
            4: ('VC', [1], [None], locations),
            5: ('DC', [1], [variant_1], locations),
            6: ('PRN', [4, 5], [None], locations)
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
        variant_1 = (0,3,'TCT','T','INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (0, 1, 'T', 'A', 'SNV', '0:T-A', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAC', [0], [None], locations),
            2: ('GCVV', [1], [None], locations),
            3: ('V', [2],[variant_1], locations),
            4: ('P', [3], [None], locations),
            5: ('VC', [2], [None], locations),
            6: ('DC', [2], [variant_2], locations),
            7: ('PRN', [5, 6], [None], locations)
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
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '', 0, 1, True)
        variant_2 = (0, 1, 'T', 'A', 'SNV', '', 0, 1, False)
        locations = locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAR', [0], [None], locations),
            2: ('GCVVV', [1], [variant_1], locations),
            3: ('P', [3], [None], locations),
            4: ('GCVVVC', [1], [None], locations),
            5: ('GCVVDC', [1], [variant_2], locations),
            6: ('PDNK', [4, 5], [None], locations),
            7: ('DCAGP', [6], [None], locations)
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
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '', 0, 1, True)
        variant_2 = (0, 1, 'T', 'A', 'SNV', '', 0, 1, False)
        locations = locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAR', [0], [None], locations),
            2: ('GCVVV', [1], [variant_1], locations),
            3: ('P', [3], [None], locations),
            4: ('GCVVVC', [1], [None], locations),
            5: ('GCVVRC', [1], [variant_2], locations),
            6: ('PDNK', [4, 5], [None], locations),
            7: ('DCAGP', [6], [None], locations)
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
        variant_1 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        variant_2 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('AER', [0], [None], locations),
            2: ('NCWH', [1], [None], locations),
            3: ('NCWV', [1], [variant_1], locations),
            4: ('STQ', [2,3], [], locations),
            5: ('Q', [4], [], locations),
            6: ('V', [4], [variant_2], locations),
            7: ('PK', [5,6], [], locations)
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
        variant_1 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        variant_2 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('AER', [0], [None], locations),
            2: ('NCWH', [1], [None], locations),
            3: ('NCWR', [1], [variant_1], locations),
            4: ('STQ', [2,3], [], locations),
            5: ('Q', [4], [], locations),
            6: ('V', [4], [variant_2], locations),
            7: ('PK', [5,6], [], locations)
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
        variant_1 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        variant_2 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('AER', [0], [None], locations),
            2: ('NCWH', [1], [None], locations),
            3: ('NCWV', [1], [variant_1], locations),
            4: ('STEEKLPAQ', [2,3], [None], locations),
            5: ('Q', [4], [None], locations),
            6: ('V', [4], [variant_2], locations),
            7: ('PK', [5,6], [None], locations)
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
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations),
            2: ('SSSSR', [1],[variant_1], locations),
            3: ('SSSIR', [1], [None], locations),
            4: ('SSXPK', [2,3], [None], locations)
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_miscleavages(self):
        """ test micleavages is handled correctly """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations),
            2: ('SSSSR', [1],[variant_1], locations),
            3: ('SSSIR', [1], [None], locations),
            4: ('SSSPK', [2,3], [None], locations),
            5: ('SSSVK', [4], [None], locations)
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR', 'SSSSRSSSPK'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_truncated(self):
        """ test tuncated peptide is handled properly """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations),
            2: ('SSSSR', [1],[variant_1], locations),
            3: ('SSSIR', [1], [None], locations),
            4: ('SSSPK', [2,3], [None], locations),
            5: ('SSSVK', [4], [None], locations)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        nodes[5].truncated = True
        peptides = graph.call_variant_peptides(2)
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR', 'SSSSKSSSSRSSSPK', 'SSSSRSSSPK'}
        self.assertEqual(received, expected)

    def test_call_peptides_check_orf(self):
        """ test calling peptides when checking for ORF """
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations, [1,None]),
            2: ('SSSSJ', [0], [None], locations, [2,None]),
            3: ('SSSSR', [1], [None], locations, [1,None]),
            4: ('SSSSR', [2], [None], locations, [2,None])
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
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations),
            2: ('SSSSR', [1],[variant_1], locations),
            3: ('SSSIR', [1], [None], locations),
            4: ('SSXPK', [2,3], [None], locations)
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        for peptide in peptides:
            self.assertTrue(peptide.description.startswith(graph.id))

    def test_call_peptides_empty_graph(self):
        """ When the graph is empty """
        data = {}
        graph, _ = create_pgraph(data, 'ENST0001')
        graph.add_stop(graph.root)
        peptides = graph.call_variant_peptides(2)
        self.assertEqual(len(peptides), 0)
