""" Test module for ThreeFrameTVG """
import copy
import unittest
from test.unit import create_variant, create_variants, \
    create_genomic_annotation, create_dna_record_dict, create_three_frame_tvg, \
    ANNOTATION_DATA, GENOME_DATA
from moPepGen import seqvar, svgraph
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation


class TestCaseThreeFrameTVG(unittest.TestCase):
    """ Test case for ThreeFrameTVG """
    def test_splice(self):
        """ Test node splice at position """
        seq = 'AAAATCCCCG'
        data = {
            1: ['AAAATCCCCG', ['RF0'], [], 0],
            2: ['AAATCCCCG', ['RF1'], [], 1],
            3: ['AATCCCCG', ['RF2'], [], 2]
        }
        graph, nodes = create_three_frame_tvg(data, seq)
        left, right = graph.splice(nodes[1], 4, 'reference')
        self.assertTrue(graph.is_out_bond_to_any_root(left))
        self.assertEqual(str(left.seq.seq), 'AAAA')
        self.assertIn(right, [e.out_node for e in left.out_edges])
        self.assertIn(left, [e.in_node for e in right.in_edges])

    def test_apply_variant_case1(self):
        """
                               C
                              / \
        AACCTTGG    ->    AACC-T-TGG
        """
        data = {
            1: ('AACCTTGG', ['RF0'], [])
        }
        graph, nodes = create_three_frame_tvg(data, 'AACCTTGG')
        variant = create_variant(4, 5, 'T', 'C', 'SNV', '')
        source, target = graph.apply_variant(nodes[1], nodes[1], variant)
        self.assertEqual(str(source.seq.seq), 'AACC')
        self.assertIs(source, target)
        seqs = {str(e.out_node.seq.seq) for e in source.out_edges}
        seqs_expected = set(['T', 'C'])
        self.assertEqual(seqs, seqs_expected)

    def test_apply_variant_after_start_codon(self):
        """ Test variant that starts right after the start codon """
        seq = 'AAAAAATGGGCCCCCTTTTT'
        data = {
            1: (seq,     ['RF0'], [], 0),
            2: (seq[1:], ['RF1'], [], 1),
            3: (seq[2:], ['RF2'], [], 2)
        }
        graph,_ = create_three_frame_tvg(data, seq)
        graph.has_known_orf = True
        graph.known_orf = [5,20]
        graph.seq.orf = FeatureLocation(start=5, end=20)
        variant = create_variant(7, 8, 'G', 'GGGCC', 'INDEL', 'INDEL-G-GGGCC')
        graph.create_variant_graph([variant], None, None, None)
        x = str(list(graph.reading_frames[0].out_edges)[0].out_node.seq.seq)
        y = 'AAAAAATGG'
        self.assertEqual(x, y)
        x = str(list(graph.reading_frames[1].out_edges)[0].out_node.seq.seq)
        y = 'AAAAATG'
        self.assertEqual(x, y)
        x = str(list(graph.reading_frames[2].out_edges)[0].out_node.seq.seq)
        y = 'AAAATG'
        self.assertEqual(x, y)

    def test_apply_fusion_case1(self):
        r""" Fusion breakpoint: acceptor 8, donor 8
             C                    A     A
            / \                  / \   / \
        AACC-T-TGGCGGTTC      TGG-G-TCC-T-TC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACC', ['RF0'], []),
            2: ('T', [1], []),
            3: ('C', [1], [(0, 'C', 'T', 'SNV', '')]),
            4: ('TGGCGGTTC', [2,3], []),
            11: ('AACC', ['RF1'], []),
            12: ('T', [11], []),
            13: ('C', [11], [(0, 'C', 'T', 'SNV', '')]),
            14: ('TGGCGGTTC', [12,13], []),
            21: ('AACC', ['RF0'], []),
            22: ('T', [21], []),
            23: ('C', [21], [(0, 'C', 'T', 'SNV', '')]),
            24: ('TGGCGGTTC', [22,23], [])
        }
        seq = 'AACCTCTGGCGGTTC'
        graph, nodes = create_three_frame_tvg(data, seq)

        attrs = {
            'ACCEPTER_GENE_ID': 'ENSG0001',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST0001.1',
            'ACCEPTER_POSITION': 20,
            'ACCEPTER_CHROM': 'chr1',
            'LEFT_INSERTION_START': None,
            'LEFT_INSERTION_END': None,
            'RIGHT_INSERTION_START': None,
            'RIGHT_INSERTION_END': None
        }
        var_fusion = create_variant(8, 9, 'C', '<FUSION>', 'Fusion',
            'FUSIONXXX', attrs=attrs)
        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (14, 15, 'G', 'A', 'SNV', '', None, 'ENST0001.1'),
            (18, 19, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        tx_variants = create_variants(var_data)
        var_data2 = {
            (25, 26, 'T', 'A', 'SNV', '' , None, 'ENSG0001')
        }
        it_variants = create_variants(var_data2)
        variant_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in tx_variants:
            variant_pool.add_transcriptional_variant(variant)
        for variant in it_variants:
            variant_pool.add_intronic_variant(variant)

        end_nodes = graph.apply_fusion([nodes[i] for i in [4,14,24]],
            var_fusion, variant_pool, genome, anno)

        for end_node in end_nodes:
            self.assertEqual(len(end_node.out_edges), 2)
            self.assertEqual(str(end_node.seq.seq), 'TGG')
            for edge in end_node.out_edges:
                if edge.type != 'reference':
                    accepter_root = edge.out_node
                    self.assertEqual(str(accepter_root.seq.seq), 'ATGG')
                    self.assertEqual(len(accepter_root.in_edges), 1)

            # also test alignment will treat the first fusioned node as reference
            graph.align_variants(end_node)
            self.assertEqual(str(end_node.seq.seq), 'TGG')

    def test_apply_fusion_intronic_donor(self):
        r""" Test fusion that the donor breakpoint is intronic
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACC', ['RF0'], []),
            2: ('T', [1], []),
            3: ('C', [1], [(0, 'C', 'T', 'SNV', '')]),
            4: ('TGGCGGTTCGACCGTG', [2,3], []),
            11: ('AACC', ['RF1'], []),
            12: ('T', [11], []),
            13: ('C', [11], [(0, 'C', 'T', 'SNV', '')]),
            14: ('TGGCGGTTCGACCGTG', [12,13], []),
            21: ('AACC', ['RF0'], []),
            22: ('T', [21], []),
            23: ('C', [21], [(0, 'C', 'T', 'SNV', '')]),
            24: ('TGGCGGTTCGACCGTG', [22,23], [])
        }
        seq = 'AACCTCTGGCGGTTC'
        graph, nodes = create_three_frame_tvg(data, seq)

        attrs = {
            'TRANSCRIPT_ID': 'ENST0001.1',
            'ACCEPTER_GENE_ID': 'ENSG0002',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST0002.1',
            'ACCEPTER_POSITION': 20,
            'ACCEPTER_CHROM': 'chr1'
        }
        var_fusion = create_variant(15, 16, 'C', '<FUSION>', 'Fusion',
            'FUSIONXXX', attrs=attrs)
        var_fusion.location.seqname= 'ENSG0001'
        var_fusion.shift_breakpoint_to_closest_exon(anno)
        var_fusion = var_fusion.to_transcript_variant(anno, genome, 'ENST0001.1')
        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0002.1'),
            (14, 15, 'G', 'A', 'SNV', '', None, 'ENST0002.1'),
            (18, 19, 'T', 'A', 'SNV', '', None, 'ENST0002.1')
        }
        tx_variants = create_variants(var_data)
        var_data2 = {
            (14, 15, 'T', 'A', 'SNV', '' , None, 'ENSG0001')
        }
        it_variants = create_variants(var_data2)
        variant_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in tx_variants:
            variant_pool.add_transcriptional_variant(variant)
        for variant in it_variants:
            variant_pool.add_intronic_variant(variant)

        end_nodes = graph.apply_fusion([nodes[i] for i in [4,14,24]],
            var_fusion, variant_pool, genome, anno)

        gene_model = anno.genes['ENSG0001']
        gene_seq = gene_model.get_gene_sequence(genome['chr1'])
        insert_seq = str(gene_seq.seq[12:15])

        for node in end_nodes:
            self.assertEqual(str(node.seq.seq), 'TG')
            received = {str(x.out_node.seq.seq) for x in node.out_edges}
            self.assertEqual(received, {insert_seq, 'GCGGTTCGACCGTG'})

    def test_apply_fusion_intronic_accepter(self):
        r""" Test fusion that the accepter breakpoint is intronic
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACC', ['RF0'], []),
            2: ('T', [1], []),
            3: ('C', [1], [(0, 'C', 'T', 'SNV', '')]),
            4: ('TGGCGGTTCGACCGTG', [2,3], []),
            11: ('AACC', ['RF1'], []),
            12: ('T', [11], []),
            13: ('C', [11], [(0, 'C', 'T', 'SNV', '')]),
            14: ('TGGCGGTTCGACCGTG', [12,13], []),
            21: ('AACC', ['RF0'], []),
            22: ('T', [21], []),
            23: ('C', [21], [(0, 'C', 'T', 'SNV', '')]),
            24: ('TGGCGGTTCGACCGTG', [22,23], [])
        }
        seq = 'AACCTCTGGCGGTTC'
        graph, nodes = create_three_frame_tvg(data, seq)

        attrs = {
            'TRANSCRIPT_ID': 'ENST0001.1',
            'ACCEPTER_GENE_ID': 'ENSG0002',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST0002.1',
            'ACCEPTER_POSITION': 14,
            'ACCEPTER_CHROM': 'chr1'
        }
        var_fusion = create_variant(20, 21, 'C', '<FUSION>', 'Fusion',
            'FUSIONXXX', attrs=attrs)
        var_fusion.location.seqname= 'ENSG0001'
        var_fusion.shift_breakpoint_to_closest_exon(anno)
        var_fusion = var_fusion.to_transcript_variant(anno, genome, 'ENST0001.1')
        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0002.1'),
            (14, 15, 'G', 'A', 'SNV', '', None, 'ENST0002.1'),
            (18, 19, 'T', 'A', 'SNV', '', None, 'ENST0002.1')
        }
        tx_variants = create_variants(var_data)
        var_data2 = {
            (16, 17, 'T', 'A', 'SNV', '' , None, 'ENSG0002')
        }
        it_variants = create_variants(var_data2)
        variant_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in tx_variants:
            variant_pool.add_transcriptional_variant(variant)
        for variant in it_variants:
            variant_pool.add_intronic_variant(variant)

        end_nodes = graph.apply_fusion([nodes[i] for i in [4,14,24]],
            var_fusion, variant_pool, genome, anno)

        gene_model = anno.genes['ENSG0002']
        gene_seq = gene_model.get_gene_sequence(genome['chr1'])
        insert_seq = str(gene_seq.seq[14:17])

        for node in end_nodes:
            self.assertEqual(str(node.seq.seq), 'TGGCG')
            received = {str(x.out_node.seq.seq) for x in node.out_edges}
            self.assertEqual(received, {insert_seq, 'GTTCGACCGTG'})

    def test_apply_fusion_intronic_both(self):
        r""" Test fusion that the both donor and accepter breakpoint are
        intronic """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ('AACC', ['RF0'], []),
            2: ('T', [1], []),
            3: ('C', [1], [(0, 'C', 'T', 'SNV', '')]),
            4: ('TGGCGGTTCGACCGTG', [2,3], []),
            11: ('AACC', ['RF1'], []),
            12: ('T', [11], []),
            13: ('C', [11], [(0, 'C', 'T', 'SNV', '')]),
            14: ('TGGCGGTTCGACCGTG', [12,13], []),
            21: ('AACC', ['RF0'], []),
            22: ('T', [21], []),
            23: ('C', [21], [(0, 'C', 'T', 'SNV', '')]),
            24: ('TGGCGGTTCGACCGTG', [22,23], [])
        }
        seq = 'AACCTCTGGCGGTTC'
        graph, nodes = create_three_frame_tvg(data, seq)

        attrs = {
            'TRANSCRIPT_ID': 'ENST0001.1',
            'ACCEPTER_GENE_ID': 'ENSG0002',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST0002.1',
            'ACCEPTER_POSITION': 14,
            'ACCEPTER_CHROM': 'chr1'
        }
        var_fusion = create_variant(15, 16, 'C', '<FUSION>', 'Fusion',
            'FUSIONXXX', attrs=attrs)
        var_fusion.location.seqname= 'ENSG0001'
        var_fusion.shift_breakpoint_to_closest_exon(anno)
        var_fusion = var_fusion.to_transcript_variant(anno, genome, 'ENST0001.1')
        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0002.1'),
            (14, 15, 'G', 'A', 'SNV', '', None, 'ENST0002.1'),
            (18, 19, 'T', 'A', 'SNV', '', None, 'ENST0002.1')
        }
        tx_variants = create_variants(var_data)
        var_data2 = {
            (16, 17, 'T', 'A', 'SNV', '' , None, 'ENSG0002')
        }
        it_variants = create_variants(var_data2)
        variant_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in tx_variants:
            variant_pool.add_transcriptional_variant(variant)
        for variant in it_variants:
            variant_pool.add_intronic_variant(variant)

        end_nodes = graph.apply_fusion([nodes[i] for i in [4,14,24]],
            var_fusion, variant_pool, genome, anno)

        gene_model = anno.genes['ENSG0001']
        gene_seq = gene_model.get_gene_sequence(genome['chr1'])
        insert_seq = str(gene_seq.seq[12:15])

        for node in end_nodes:
            self.assertEqual(str(node.seq.seq), 'TG')
            received = {str(x.out_node.seq.seq) for x in node.out_edges}
            self.assertEqual(received, {insert_seq, 'GCGGTTCGACCGTG'})

    def test_apply_variant_inframe(self):
        """ apply_variant """
        data = {
            1: ['AAAAAT', ['RF0'], [], 0],
            2: ['AAAAT', ['RF1'], [], 1],
            3: ['AAAT', ['RF2'], [], 2]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAAAAT')
        variant = create_variant(3, 4, 'A', 'T', 'SNV', '')
        source, target = graph.apply_variant(nodes[1], nodes[1], variant)
        self.assertIs(source, target)
        out_nodes = {str(e.out_node.seq.seq) for e in source.out_edges}
        expected = {'A', 'T'}
        self.assertEqual(out_nodes, expected)

    def test_apply_variant_inframe_case2(self):
        """ apply inframe variant, where the reference node start position
        is the same as the variant position. """
        var_data = [
            [4, 5, 'T', 'A', 'SNV', ''],
            [5, 6, 'C', 'G', 'SNV', '']
        ]
        variants = create_variants(var_data)
        data = {
            1: ['AAAA', ['RF0'], [], 0],
            2: ['T', [1], [], 0],
            3: ['A', [1], [(0,'T','A','SNV','')], 0],
            4: ['CTTAC', [2,3], [], 0]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAAAAT')
        loc = MatchedLocation(
            query=FeatureLocation(0,5),
            ref=FeatureLocation(5,10)
        )
        nodes[4].seq.locations = [loc]
        source, target = graph.apply_variant(nodes[4], nodes[4], variants[1])
        self.assertIs(source, target)
        self.assertTrue(nodes[2].is_inbond_of(source))
        self.assertTrue(nodes[3].is_inbond_of(source))
        self.assertEqual(str(source.seq.seq), 'C')
        new_node = [x for x in nodes[2].out_edges if x.out_node is not source][0]\
            .out_node
        self.assertEqual(new_node.seq.seq, 'G')

    def test_apply_variant_deletion(self):
        """ apply_variant """
        data = {
            1: ['AAAAAT', ['RF0'], [], 0],
            2: ['AAAAT', ['RF1'], [], 1],
            3: ['AAAT', ['RF2'], [], 2]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAAAAT')
        variant = create_variant(3, 5, 'AA', 'A', 'INDEL', '')
        source, target = graph.apply_variant(nodes[1], nodes[2], variant)
        self.assertTrue(graph.is_out_bond_to_any_root(source))
        self.assertTrue(graph.is_out_bond_to_any_root(source))
        self.assertEqual(str(source.seq.seq), 'AAA')
        self.assertEqual(str(target.seq.seq), 'AAAA')
        out_nodes = {str(e.out_node.seq.seq) for e in source.out_edges}
        self.assertEqual(out_nodes, {'A', 'AAT'})

    def test_create_variant_graph_in_frame(self):
        """ apply_variant """
        data = {
            1: ['AAATAAATAAAT', ['RF0'], [], 0],
            2: ['AATAAATAAAT', ['RF1'], [], 1],
            3: ['ATAAATAAAT', ['RF2'], [], 2]
        }

        var_data = [
            (3, 4, 'T', 'A', 'SNV', ''),
            (6, 7, 'A', 'T', 'SNV', '')
        ]

        graph, _ = create_three_frame_tvg(data, 'AAATAAATAAAT')
        variants = create_variants(var_data)
        graph.create_variant_graph(variants, None, None, None)
        received = {str(list(n.out_edges)[0].out_node.seq.seq)
            for n in graph.reading_frames}
        expected = {'AAA', 'AA', 'A'}
        self.assertEqual(received, expected)
        received = {str(e.out_node.seq.seq) for e in \
            list(graph.reading_frames[0].out_edges)[0].out_node.out_edges}
        expected = {'A', 'T'}
        self.assertEqual(received, expected)

    def test_create_variant_graph_frameshifting(self):
        """ apply_variant """
        data = {
            1: ['AAATAAATAAAT', ['RF0'], [], 0],
            2: ['AATAAATAAAT',  ['RF1'], [], 1],
            3: ['ATAAATAAAT',   ['RF2'], [], 2]
        }

        var_data = [
            (3, 5, 'TA', 'T', 'INDEL', ''),
            (7, 8, 'T', 'TC', 'INDEL', '')
        ]

        graph, _ = create_three_frame_tvg(data, 'AAATAAATAAAT')
        variants = create_variants(var_data)
        graph.create_variant_graph(variants, None, None, None)
        received = {str(list(n.out_edges)[0].out_node.seq.seq) \
            for n in graph.reading_frames}
        expected = {'AAA', 'AA', 'A'}
        self.assertEqual(received, expected)

    def test_apply_insertion_case1(self):
        """ apply_insertion """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ['AAATAAATAAAT', ['RF0'], [], 0],
            2: ['AATAAATAAAT',  ['RF1'], [], 1],
            3: ['ATAAATAAAT',   ['RF2'], [], 2]
        }
        ins_attrs = {
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        ins_var = create_variant(3,4,'T','<INS>','Insertion','',ins_attrs)

        var_data = [(10, 11, 'A', 'T', 'SNV', '')]
        variants = create_variants(var_data)
        var_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in variants:
            var_pool.add_transcriptional_variant(variant)

        graph, nodes = create_three_frame_tvg(data, 'AAATAAATAAAT')
        cursors = [nodes[1], nodes[2], nodes[3]]
        graph.apply_insertion(cursors, ins_var, var_pool, genome, anno)

        received = {str(list(n.out_edges)[0].out_node.seq.seq)
            for n in graph.reading_frames}
        self.assertEqual(received, {'A', 'AA', 'AAA'})

    def test_apply_insertion_case2(self):
        r"""
                                       T    A
                                      / \  / \
                                    CC-C-TA-T-G
                                   /           \
        AACCTGGCGGTTC     ->    AAC-C-----------TGGCGGTTC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
             1: ('AACCTGGCGGTTC', ['RF0'], [], 0),
            11: ( 'ACCTGGCGGTTC', ['RF1'], [], 1),
            21: (  'CCTGGCGGTTC', ['RF2'], [], 2)
        }
        seq = 'AACCTGGCGGTTC'
        graph, nodes = create_three_frame_tvg(data, seq)

        attrs = {
            'GENE_ID': 'ENSG0001',
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        var_insertion = create_variant(3,4,'C','<INS>','Insertion','', attrs)

        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (8, 9, 'C', 'T', 'SNV', '', None, 'ENST0001.1'),
            (11, 12, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }
        variants = create_variants(var_data)
        variant_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in variants:
            variant_pool.add_transcriptional_variant(variant)

        cursors = [nodes[i] for i in [1,11,21]]
        end_nodes = graph.apply_insertion(cursors, var_insertion,
            variant_pool, genome, anno)

        def filter_func(x):
            return lambda xx:str(xx.out_node.seq.seq) == x

        received = {str(node.seq.seq) for node in end_nodes}
        expected = {'AAC', 'AC', 'C'}
        self.assertEqual(received, expected)

        for node in end_nodes:
            node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
            self.assertEqual(node_seqs, {'C', 'CC'})

            node = list(filter(filter_func('CC'), node.out_edges))[0].out_node
            node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
            self.assertEqual(node_seqs, {'C', 'T'})

            node = list(filter(filter_func('C'), node.out_edges))[0].out_node
            node = list(node.out_edges)[0].out_node
            node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
            self.assertEqual(node_seqs, {'T', 'A'})

    def test_apply_insertion_case3(self):
        """ apply_insertion with end_inclusion """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            1: ['AAATAAATAAAT', ['RF0'], [], 0],
            2: ['AATAAATAAAT',  ['RF1'], [], 1],
            3: ['ATAAATAAAT',   ['RF2'], [], 2]
        }
        gene_id = 'ENSG0001'
        ins_attrs = {
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        ins_var = create_variant(2,3,'T','<INS>','Insertion','',ins_attrs)
        gene_seq = anno.genes[gene_id].get_gene_sequence(genome['chr1'])
        ins_var.to_end_inclusion(gene_seq)

        var_pool = seqvar.VariantRecordPool(anno=anno)
        graph, nodes = create_three_frame_tvg(data, 'AAATAAATAAAT')
        cursors = [nodes[1], nodes[2], nodes[3]]
        graph.apply_insertion(cursors, ins_var, var_pool, genome, anno)

        received = {str(list(n.out_edges)[0].out_node.seq.seq)
            for n in graph.reading_frames}
        self.assertEqual(received, {'A', 'AA', 'AAA'})

        for root in graph.reading_frames:
            node = root.get_out_nodes()[0]
            for x in node.get_out_nodes():
                if x.variants:
                    self.assertEqual(x.seq.seq, 'CCTATGT')
                else:
                    self.assertEqual(x.seq.seq, 'T')

    def test_apply_substitution_case1(self):
        r"""                         T
                                    / \
                                   C-C-TATG
                                  /         \
        AACCTGGCGGTTC    ->    AAC-CTGGC-----GGTTC
        """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
             1: ('AACCTGGCGGTTC', ['RF0'], [], 0),
            11: ( 'ACCTGGCGGTTC', ['RF1'], [], 1),
            21: (  'CCTGGCGGTTC', ['RF2'], [], 2)
        }
        seq = 'AACCTGGCGGTTC'
        graph, nodes = create_three_frame_tvg(data, seq)

        attrs = {
            'GENE_ID': 'ENSG0001',
            'START': 3,
            'END': 8,
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        var_sub = create_variant(3,8,'C','<SUB>','Substitution','', attrs)

        var_data = {
            ( 3,  4, 'T', 'A', 'SNV', '', None, 'ENST0001.1'),
            (8, 9, 'C', 'T', 'SNV', '', None, 'ENST0001.1'),
            (19, 20, 'T', 'A', 'SNV', '', None, 'ENST0001.1')
        }

        variants = create_variants(var_data)
        variant_pool = seqvar.VariantRecordPool(anno=anno)
        for variant in variants:
            variant_pool.add_transcriptional_variant(variant)

        cursors = [nodes[i] for i in [1,11,21]]
        end_nodes = graph.apply_substitution(cursors, var_sub, variant_pool,
            genome, anno)

        received = {str(node.seq.seq) for node in end_nodes}
        expected = {'AAC' ,'AC', 'C'}
        self.assertEqual(received, expected)

        for node in end_nodes:
            node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
            self.assertEqual(node_seqs, {'CTGGC', 'C'})

            func = lambda x:str(x.out_node.seq.seq) == 'C'
            node = list(filter(func, node.out_edges))[0].out_node
            node_seqs = {str(edge.out_node.seq.seq) for edge in node.out_edges}
            self.assertEqual(node_seqs, {'C', 'T'})

    def test_apply_deletion_case1(self):
        r"""
                                   C-------
                                  /        \
        AACCTGGCGGTTC    ->    AAC-CTGGC----GGTTC
        """
        seq = 'AACCTGGCGGTTC'
        data = {
             1: [seq,     ['RF0'], [], 0],
            11: [seq[1:], ['RF0'], [], 1],
            21: [seq[2:], ['RF0'], [], 2]
        }
        del_var = create_variant(3,8,'C','<DEL>','Deletion','')
        del_var.attrs['GENE_ID'] = 'ENSG0001'
        graph, nodes = create_three_frame_tvg(data, seq)

        cursors = [nodes[i] for i in [1,11,21]]

        frames_shifted = (-(8 - 3)) % 3
        for i in range(3):
            j = (i + frames_shifted) % 3
            cursors[i], cursors[j] = graph.apply_variant(cursors[i],
                cursors[j], del_var)

        received = {str(node.seq.seq) for node in cursors}
        expected = {'AAC' ,'AC', 'C'}
        self.assertEqual(received, expected)

        for node in cursors:
            self.assertEqual(len(node.out_edges), 2)
            for edge in node.out_edges:
                x = edge.out_node
                if edge.type == 'reference':
                    self.assertEqual(str(x.seq.seq), 'CTGGC')
                else:
                    self.assertEqual(str(x.seq.seq), 'C')
                    self.assertEqual(len(x.in_edges), 1)
                    self.assertEqual(len(x.out_edges), 1)

    def test_truncate_left_variant_right(self):
        """ test truncate left when variant on the right"""
        data = {
            1: ('ATGG',  ['RF0'], []),
            2: ('CCATG', [1],     [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1],     []),
            4: ('CCCT',  [2,3],   [])
        }
        seq = 'ATGGCCCTGCCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        left = nodes[2].truncate_left(2)
        self.assertEqual(str(left.seq.seq),'CC')
        self.assertEqual(len(left.variants), 0)
        self.assertEqual(str(nodes[2].seq.seq),'ATG')
        self.assertEqual(len(nodes[2].variants), 1)
        self.assertEqual(nodes[2].variants[0].location.start, 0)
        self.assertTrue(nodes[2].is_inbond_of(nodes[4]))
        self.assertFalse(left.in_edges)
        self.assertFalse(left.out_edges)

    def test_truncate_left_variant_left(self):
        """ test truncate left when variant on the left """
        data = {
            1: ('ATGG',  ['RF0'], []),
            2: ('CCATG', [1],     [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1],     []),
            4: ('CCCT',  [2,3],   [])
        }
        seq = 'ATGGCCCTGCCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        left = nodes[2].truncate_left(3)
        self.assertEqual(str(left.seq.seq),'CCA')
        self.assertEqual(len(left.variants), 1)
        self.assertEqual(left.variants[0].location.start, 2)
        self.assertEqual(str(nodes[2].seq.seq),'TG')
        self.assertEqual(len(nodes[2].variants), 0)

    def test_truncate_right_variant_left(self):
        """ test truncate left when variant on the left """
        data = {
            1: ('ATGG',  ['RF0'], []),
            2: ('CCATG', [1],     [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1],     []),
            4: ('CCCT',  [2,3],   [])
        }
        seq = 'ATGGCCCTGCCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        right = nodes[2].truncate_right(3)
        self.assertEqual(str(nodes[2].seq.seq),'CCA')
        self.assertEqual(len(nodes[2].variants), 1)
        self.assertEqual(nodes[2].variants[0].location.start, 2)
        self.assertEqual(str(right.seq.seq),'TG')
        self.assertEqual(len(right.variants), 0)
        self.assertTrue(nodes[1].is_inbond_of(nodes[2]))
        self.assertFalse(right.in_edges)
        self.assertFalse(right.out_edges)

    def test_truncate_right_variant_right(self):
        """ test truncate left when variant on the right """
        data = {
            1: ('ATGG',  ['RF0'], []),
            2: ('CCATG', [1],     [(2, 'G', 'A', 'SNV', '')]),
            3: ('CCCTG', [1],     []),
            4: ('CCCT',  [2,3],   [])
        }
        seq = 'ATGGCCCTGCCCT'
        _, nodes = create_three_frame_tvg(data, seq)
        right = nodes[2].truncate_right(2)
        self.assertEqual(str(nodes[2].seq.seq),'CC')
        self.assertEqual(len(nodes[2].variants), 0)
        self.assertEqual(str(right.seq.seq),'ATG')
        self.assertEqual(len(right.variants), 1)
        self.assertEqual(right.variants[0].location.start, 0)

    def test_expand_alignments_case1(self):
        """ Variant alignments are expanded """
        seq = 'ATGGTCTGCCCTCTGAACTGA'
        data = {
            1: ('ATG', ['RF0'], []),
            2: ('G', [1], []),
            3: ('T', [1], [(3, 'G', 'T', 'SNV', '3:G-T')]),
            4: ('A', [1], []),
            5: ('TCTGCCCTCTGAACTGA', [2,3,4], [])
        }
        graph, nodes = create_three_frame_tvg(data, seq)
        node = nodes[1]
        graph.expand_alignments(node)
        variant_site_nodes = [edge.out_node for edge in node.out_edges]
        variant_site_seqs = [str(node.seq.seq) for node in variant_site_nodes]
        self.assertEqual(set(['TTC', 'ATC', 'GTC']), set(variant_site_seqs))

        seqs = [edge.in_node.seq.seq for edge in node\
            .get_reference_next().get_reference_next().in_edges]
        self.assertEqual(set(['TTC', 'ATC', 'GTC']), set(seqs))

    def test_merge_with_outbounds(self):
        """
            ATGG-T-G-CCCT
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T',    [1],     []),
            3: ('G',    [2],     []),
            4: ('CCCT', [3],     [])
        }
        seq = 'ATGGTGCCCT'
        graph, nodes = create_three_frame_tvg(data, seq)
        returned_nodes = graph.merge_with_outbonds(nodes[3])
        self.assertEqual(str(nodes[4].seq.seq), 'GCCCT')
        self.assertTrue(nodes[2].is_inbond_of(nodes[4]))
        self.assertFalse(nodes[2].is_inbond_of(nodes[3]))
        self.assertFalse(nodes[3].is_inbond_of(nodes[4]))
        self.assertEqual(returned_nodes[0], nodes[4])

    def test_expand_alignments_case2(self):
        r""" Tests that the node after variants are truncated.
                 A
                / \
            ATGG-G-CCCT
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('A',    [1],     [(0, 'G', 'A', 'SNV', '')]),
            3: ('G',    [1],     []),
            4: ('CCCT', [2,3],   [])
        }
        seq = 'ATGGAGCCCT'
        graph, nodes = create_three_frame_tvg(data, seq)
        end_nodes = graph.expand_alignments(nodes[1])
        self.assertEqual(str(nodes[4].seq.seq), 'CCT')
        self.assertIs(nodes[4], end_nodes[0])


    def test_align_variants(self):
        """ find_known_orf wihtout mutation """
        data = {
            1: ['AAAT', ['RF0'], []],
            2: ['A', [1], []],
            3: ['G', [1], [(0, 'A', 'G', 'SNV', '')]],
            4: ['G', [2,3], []],
            5: ['T', [2], [(0, 'G', 'T', 'SNV', '')]],
            6: ['C', [4,5], []]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAATAAATAAAT')
        graph.align_variants(nodes[1])
        self.assertEqual(len(nodes[1].out_edges), 3)
        received = {str(e.out_node.seq.seq) for e in nodes[1].out_edges}
        expected = {'GG', 'AG', 'AT'}
        self.assertEqual(received, expected)

    def test_align_variants_bridge(self):
        r""" find_known_orf wihtout mutation

        0: AAA-TA-G            0: AAA-TA-G
              \                      \
               T                      TG
                \                       \
               G |      ->            GG |
              / \|                   /  \|
        1: AAT-A-G-C           1: AAT-AG-C
                \ /                  \  /
                 T                    AT
        """
        data = {
            1:  ['AAA', ['RF0'], []],
            2:  ['TA', [1], []],
            3:  ['T', [1], [(0, 'TA', 'T', 'INDEL', '')]],
            4:  ['G', [2], []],
            11: ['AAT', ['RF1'], [], 1],
            12: ['A', [11], []],
            13: ['G', [11], [(0, 'A', 'G', 'SNV', '')]],
            14: ['G', [3, 12,13], []],
            15: ['T', [12], [(0, 'G', 'T', 'SNV', '')]],
            16: ['C', [14,15], []]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAATAAATAAAT')
        graph.align_variants(nodes[11])
        self.assertEqual(len(nodes[11].out_edges), 3)
        received = {str(e.out_node.seq.seq) for e in nodes[11].out_edges}
        expected = {'GG', 'AG', 'AT'}
        self.assertEqual(received, expected)
        self.assertEqual(len(nodes[16].in_edges), 4)
        received = {str(e.in_node.seq.seq) for e in nodes[16].in_edges}
        expected = {'GG', 'AG', 'AT', 'TG'}
        self.assertEqual(received, expected)

    def test_align_variant_sub(self):
        r""" when there is a subgraph out node

        0: AAAAA-TA-GTT            0: AAAAA-TA-GTT
                \                          \
                 TGCTGC-T-GC                TGCTGC-T-GC
                       \ /                        \ /
                        A                          A
        """
        data = {
            1:  ['AAAAA', ['RF0'], []],
            2:  ['TA', [1], []],
            3:  ['GTT', [2], []]
        }
        graph, nodes = create_three_frame_tvg(data, 'AAAAATAGTT', 'ENST01')
        data = {
            11: ['TGCTGC', ['RF0'], []],
            12: ['T', [11], []],
            13: ['A', [12], [(0, 'A', 'G', 'SNV', '')]],
            14: ['GC', [12, 13], []]
        }
        graph2, nodes2 = create_three_frame_tvg(data, 'TGCTGCTAGC')
        for edge in copy.copy(nodes2[11].in_edges):
            graph2.remove_edge(edge)
        graph2.add_edge(nodes[1], nodes2[11], 'variant_start')
        graph.align_variants(nodes[1])
        out_node_seqs = {x.out_node.seq.seq for x in nodes[1].out_edges}
        expected = {'TA', 'TGCTGC'}
        self.assertEqual(out_node_seqs, expected)

    def test_align_variantes_max_variants(self):
        r""" Tests for aligning variants in heavily mutated local region

            G
           / \
          | T |   T   A
          |/ \|  / \ / \
        AAT-A-C-C-C-T-T-G
           \ /
            C

        """
        var_data = [
            (0, 'A', 'T', 'SNV', ''),
            (0, 'A', 'G', 'SNV', ''),
            (0, 'A', 'C', 'SNV', ''),
            (0, 'C', 'G', 'SNV', ''),
            (0, 'C', 'G', 'SNV', ''),
            (0, 'T', 'A', 'SNV', ''),
        ]
        data = {
            1:  ['AAT', ['RF0'], [], 1],
            2:  ['A', [1], []],
            3:  ['T', [1], [var_data[0]]],
            4:  ['G', [1], [var_data[1]]],
            5:  ['C', [1], [var_data[2]]],
            6:  ['C', [2,3,4,5], []],
            7:  ['C', [6], []],
            8:  ['C', [6], [var_data[3]]],
            9:  ['T', [7,8], []],
            10: ['T', [9], []],
            11: ['A', [9], [var_data[4]]],
            12: ['G', [10,11], []]
        }

        graph, nodes = create_three_frame_tvg(data, 'AATACCTTG')
        graph.cleavage_params.max_variants_per_node = 2
        graph.align_variants(nodes[1])
        for edge in nodes[1].out_edges:
            out_node = edge.out_node
            self.assertTrue(len(out_node.variants) <=
                graph.cleavage_params.max_variants_per_node)

    def test_expand_alignments(self):
        r""" find_known_orf wihtout mutation

             GG                  TGG
            /  \                /   \
        AAAT-AG-CCC    ->    AAA-TAG-CCC
            \  /                \   /
             AT                  TAT
        """
        data = {
            1: ['AAAT', ['RF0'], []],
            2: ['AG', [1], []],
            3: ['GG', [1], [(0, 'A', 'G', 'SNV', '')]],
            4: ['AT', [1], [(1, 'G', 'T', 'SNV', '')]],
            5: ['CCC', [2,3,4], []]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAATAAATAAAT')
        graph.expand_alignments(nodes[1])
        self.assertEqual(str(nodes[1].seq.seq), 'AAA')
        received = {str(e.out_node.seq.seq) for e in nodes[1].out_edges}
        expected = {'TGG', 'TAG', 'TAT'}
        self.assertEqual(received, expected)

    def test_expand_alignment_bridge(self):
        r""" find_known_orf wihtout mutation

            0: AAA-TA-G           0: AAA-TA-G
                  \                     \
                   TG-                   TG-
                      \                     \
                    GG |     ->          TGG |
                   /  \|                /   \|
            1: AAAT-AG-CCC        1: AAA-TAG-CCC
                   \  /                 \   /
                    AT                   TAT
        """
        data = {
            1:  ['AAA', ['RF0'], []],
            2:  ['TA', [1], []],
            3:  ['TG', [1], [(0, 'TA', 'T', 'INDEL', '')]],
            4:  ['G', [2], []],
            11: ['AAAT', ['RF1'], []],
            12: ['AG', [11], []],
            13: ['GG', [11], [(0, 'A', 'G', 'SNV', '')]],
            14: ['AT', [11], [(1, 'G', 'T', 'SNV', '')]],
            15: ['CCC', [12,13,14, 3], []]
        }

        graph, nodes = create_three_frame_tvg(data, 'AAATAAATAAAT')
        graph.expand_alignments(nodes[11])
        self.assertEqual(len(nodes[11].out_edges), 3)
        self.assertEqual(str(nodes[11].seq.seq), 'AAA')
        received = {str(e.out_node.seq.seq) for e in nodes[11].out_edges}
        expected = {'TGG', 'TAG', 'TAT'}
        self.assertEqual(received, expected)
        received = {str(e.in_node.seq.seq) for e in nodes[15].in_edges}
        expected = {'TG', 'TGG', 'TAG', 'TAT'}
        self.assertEqual(received, expected)

    def test_expand_variants_sub(self):
        r""" when there is a subgraph out node

        0: AAAAA-TA-GTT            0: AAA-AATAGT-T
                \                        \
                 TGCTGC-T-GC              AATGCTGC-T-GC
                       \ /                        \ /
                        A                          A
        """
        data = {
            1:  ['AAAAA', ['RF0'], []],
            2:  ['TA', [1], []],
            3:  ['GTT', [2], []]
        }
        graph, nodes = create_three_frame_tvg(data, 'AAAAATAGTT', 'ENST01')
        data = {
            11: ['TGCTGC', ['RF0'], []],
            12: ['T', [11], []],
            13: ['A', [12], [(0, 'A', 'G', 'SNV', '')]],
            14: ['GC', [12, 13], []]
        }
        graph2, nodes2 = create_three_frame_tvg(data, 'TGCTGCTAGC')
        for edge in copy.copy(nodes2[11].in_edges):
            graph2.remove_edge(edge)
        graph2.add_edge(nodes[1], nodes2[11], 'variant_start')
        graph.expand_alignments(nodes[1])
        self.assertEqual(nodes[1].seq.seq, 'AAA')
        out_node_seqs = {x.out_node.seq.seq for x in nodes[1].out_edges}
        expected = {'AATAGT', 'AATGCTGC'}
        self.assertEqual(out_node_seqs, expected)

    def test_translate_mrna_end_nf(self):
        """ Test that the UTR region is not translated with mrna_end_NF """
        data = {
            1: ['ATGAAA', ['RF0'], []],
            2: ['AAAAAA', [1], []],
            3: ['AAAAAA', [2], []]
        }
        graph, _ = create_three_frame_tvg(data, 'ATGAAAAAAAAAAAAAAA')
        graph.mrna_end_nf = True
        graph.has_known_orf = True
        graph.seq.orf = FeatureLocation(start=0, end=15)
        pgraph = graph.translate()
        x = [x for x in pgraph.root.out_nodes if x.seq.seq == 'MK'][0]
        node = list(list(x.out_nodes)[0].out_nodes)[0]
        self.assertEqual(str(node.seq.seq), 'K')

    def test_fusion_breakpoint_end_of_transcript(self):
        """ Test case for fusion that the donor breakpoint is the end of the
        transcript """
        genome = create_dna_record_dict(GENOME_DATA)
        anno = create_genomic_annotation(ANNOTATION_DATA)

        fusion_attrs = {
            'GENE_ID': 'ENSG0001',
            'GENE_SYMBOL': 'SYMBO1',
            'GENOMIC_POSITION': 'chr1-36:36',
            'ACCEPTER_GENE_ID': 'ENSG0002',
            'ACCEPTER_TRANSCRIPT_ID': 'ENST0002.1',
            'ACCEPTER_SYMBOL': 'SYMB2',
            'ACCEPTER_POSITION': 17,
            'ACCEPTER_GENOMIC_POSITION': 'chr1-78:78',
            'LEFT_INSERTION_START': None,
            'LEFT_INSERTION_END': None,
            'RIGHT_INSERTION_START': None,
            'RIGHT_INSERTION_END': None
        }
        fusion = create_variant(21, 22, 'A', '<FUS>', 'Fusion', 'FUSION-XXX',
            fusion_attrs, 'ENST0001.1')
        variant_pool = seqvar.VariantRecordPoolOnDisk(
            pointers=None, anno=anno, genome=genome
        )
        tx_model = anno.transcripts['ENST0001.1']
        tx_seq = tx_model.get_transcript_sequence(genome['chr1'])
        tgraph = svgraph.ThreeFrameTVG(tx_seq, 'ENST0001.1')
        tgraph.init_three_frames()
        tgraph.create_variant_graph([fusion], variant_pool, genome, anno)
        node = list(tgraph.reading_frames[0].out_edges)[0].out_node
        self.assertTrue(node.get_out_nodes()[0].level == 1)
        tgraph.fit_into_codons()
        pgraph = tgraph.translate()
        peptides = pgraph.call_variant_peptides()
        expected = {'MGPSFCEF', 'GPSFCEF'}
        self.assertEqual({str(x.seq) for x in peptides}, expected)

    def test_find_farthest_node_with_overlap_case1(self):
        r"""
                 T--
                /   \
            ATGG-TCT-G-CCCT
                    \ /
                     A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], [])
        }
        seq = 'ATGGTCTGACCCT'
        graph, nodes = create_three_frame_tvg(data, seq)
        node = graph.find_farthest_node_with_overlap(nodes[1])
        self.assertIs(node, nodes[6])

    def test_find_farthest_node_with_overlap_case2(self):
        r""" Test case for the last node is too short. Should include the
        following node.
                 T--   T
                /   \ / \
            ATGG-TCT-G-C-C-CT
                    \ /
                     A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCT', 'T', 'INDEL', '')]),
            3: ('TCT', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('T', [4], [(0, 'C', 'T', 'SNV', '')]),
            7: ('C', [4,5], []),
            8: ('C', [6,7], []),
            9: ('CT', [8], [])
        }
        seq = 'ATGGTCTGCCCT'
        graph, nodes = create_three_frame_tvg(data, seq)
        node = graph.find_farthest_node_with_overlap(nodes[1])
        print(node.seq.seq)
        self.assertIs(node, nodes[9])

    def test_find_farthest_node_with_overlap_case3(self):
        r""" For a bubble that do not merge, should return None. This would  be
        the case of mutation at the last nucleotide.
                 C
                /
            ATGG-T
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], []),
            3: ['C', [1], [(0, 'T', 'C', 'SNV', '')]]
        }
        seq = 'ATGGT'
        graph, nodes = create_three_frame_tvg(data, seq)
        node = graph.find_farthest_node_with_overlap(nodes[1])
        self.assertEqual(str(node.seq.seq), 'T')

    def test_find_farthest_node_with_overlap_case5_with_branch(self):
        r"""
                 T-G-CCCT
                /
            ATGG-TCTGAC-G-CCCT
                       \ /
                        A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            5: ('TCTGAC', [1], []),
            6: ('G', [5], []),
            7: ('A', [5], [(0, 'G', 'T', 'SNV', '')]),
            8: ('CCCT', [6, 7], [])
        }
        seq = 'ATGGTCTGACGCCCT'
        graph, nodes = create_three_frame_tvg(data, seq)
        node = graph.find_farthest_node_with_overlap(nodes[1])
        self.assertIs(node, nodes[5])

    def test_find_farthest_node_with_exclusive_outbond(self):
        r"""
                 T--
                /    \
            ATGG-TCTC-G-CCCT-GTTGGCCC
                     \ /
                      A
        """
        data = {
            1: ('ATGG', ['RF0'], []),
            2: ('T', [1], [(0, 'TCTGAC', 'T', 'INDEL', '')], True),
            3: ('TCTC', [1], []),
            4: ('G', [2,3], []),
            5: ('A', [3], [(0, 'G', 'T', 'SNV', '')]),
            6: ('CCCT', [4,5], []),
            7: ('GTTGGCCC', [5,6], []),
        }
        seq = 'ATGGTCTCGCCCTGTTGGCCC'
        graph, nodes = create_three_frame_tvg(data, seq)
        node = graph.find_farthest_node_with_overlap(nodes[1])
        self.assertIs(node, nodes[7])
