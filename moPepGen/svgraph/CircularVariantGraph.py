""" Module for CircularVariantGraph, a directed cyclic graph, for circRNA etc.
"""
from typing import Union, List, Deque, Tuple, Dict
import copy
import re
from collections import deque
from moPepGen import svgraph, dna, seqvar


class CircularVariantGraph(svgraph.TranscriptVariantGraph):
    """ Defines a directed cyclic graph for circular nucleotide molecules such
    ass circRNA and the variants associated with it.

    Attributes:
        seq (DNASeqRecordWithCoordinates): The original sequence of the
                transcript (reference).
        transcript_id (str): The transcript ID that the circRNA is associated
            with
        attrs (dict): additional attributes
    """
    def __init__(self, seq:Union[dna.DNASeqRecordWithCoordinates,None],
            transcript_id:str, attrs:dict=None):
        """ Construct a CircularVariantGraph

        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the
                    transcript (reference).
            transcript_id (str): The transcript ID that the circRNA is
                associated with.
            attrs (dict): additional attributes

        """
        self.seq = seq
        if self.seq and not self.seq.locations:
            self.add_default_sequence_locations()
        node = svgraph.TVGNode(seq)
        self.root = node
        # self.add_null_root()
        self.add_edge(node, self.root, _type='reference')
        self.transcript_id = transcript_id
        self.attrs = attrs

    def create_variant_graph(self, variants: List[seqvar.VariantRecord]):
        """ Apply a list of variants to the graph. Variants not in the
        range are ignored. Variants at the first nucleotide of each fragment
        of the sequence are also ignored, because it causes the exon splice
        site to be changed.
        """
        filtered_variants = []
        for variant in variants:
            for location in self.seq.locations:
                if variant.location.start > location.ref.start and \
                    variant.location.end < location.ref.end:
                    filtered_variants.append(variant)
                    break
        return super().create_variant_graph(filtered_variants)

    def create_branch(self, node:svgraph.TVGNode) -> svgraph.TVGNode:
        """ Create a branch of a given node and all its downstream nodes. The
        returned node is a leading node of a directed acyclic graph, and can be
        added into a TranscriptVariantGraph.

        Args:
            node (svgraph.TVGNode): The node to create a branch
        """
        new_node = svgraph.TVGNode(
            seq=node.seq,
            variants=copy.copy(node.variants),
            frameshifts=copy.copy(node.frameshifts)
        )

        T = Deque[Tuple[svgraph.TVGNode, svgraph.TVGNode, int]]
        queue:T = deque([(node, new_node, 0)])
        visited:Dict[Tuple(svgraph.TVGNode, int), svgraph.TVGNode] = {}

        while queue:
            source, target, _round = queue.pop()
            if source is node:
                _round += 1
            if _round >= 4:
                break
            if (source, _round) in visited:
                continue
            for edge in source.out_edges:
                source_out_node = edge.out_node
                if (source_out_node, _round) in visited:
                    new_out_node = visited[(source_out_node, _round)]
                else:
                    frameshifts = copy.copy(source_out_node.frameshifts)
                    frameshifts.update(source.frameshifts)
                    new_out_node = svgraph.TVGNode(
                        seq=source_out_node.seq,
                        variants=copy.copy(source_out_node.variants),
                        frameshifts=frameshifts
                    )
                    visited[(source_out_node, _round)] = new_out_node
                new_edge = svgraph.TVGEdge(target, new_out_node,
                    _type=edge.type)
                target.out_edges.add(new_edge)
                new_out_node.in_edges.add(new_edge)
                queue.appendleft((edge.out_node, new_out_node, _round))
        return new_node

    def find_orf_in_outbound_nodes(self, node:svgraph.TVGNode,
            carry_over:dna.DNASeqRecordWithCoordinates
            ) -> Tuple[svgraph.TVGNode, List[svgraph.TVGNode]]:
        """ Look for potential start codon int the out bound nodes. """
        branches = []
        downstream = node.find_farthest_node_with_overlap()

        for edge in node.out_edges:
            seq = carry_over + edge.out_node.seq + downstream.seq[:2]
            indices = [m.start() for m in re.finditer('ATG', str(seq.seq))]
            for index in indices:
                branch = self.create_branch(edge.out_node)
                branch.seq = (carry_over + edge.out_node.seq)[index:]
                filtered_variants = []
                for variant in branch.variants:
                    if int(variant.location.start) + len(carry_over.seq) - index > 0:
                        filtered_variants.append(variant)
                branch.variants = filtered_variants
                branch.branch = edge.type != 'reference'
                branches.append(branch)
        return downstream, branches


    def find_all_orfs(self) -> svgraph.TranscriptVariantGraph:
        """ Find all ORFs """
        tvg = svgraph.TranscriptVariantGraph(None, self.attrs['id'])

        self.align_all_variants()

        queue:Deque[svgraph.TVGNode] = deque([self.root])

        while queue:
            cur = queue.pop()

            indices = [m.start() for m in re.finditer('ATG', str(cur.seq.seq))]
            for index in indices:
                branch = self.create_branch(cur)
                branch.seq = branch.seq[index:]
                _type = 'variant_start' if branch.branch else 'reference'
                tvg.add_edge(tvg.root, branch, _type)

            index = indices[-1] if indices else -2
            carry_over = cur.seq[index:]

            main, branches = self.find_orf_in_outbound_nodes(cur, carry_over)
            for branch in branches:
                _type = 'variant_start' if branch.branch else 'reference'
                tvg.add_edge(tvg.root, branch, _type)

            if main.is_inbond_of(self.root):
                queue.append(main)
        return tvg


    def align_all_variants(self) -> None:
        """ Align all variants """
        node = self.root
        while True:
            if len(node.out_edges) > 1:
                self.align_variants(node, branch_out_frameshifting=False)
                node = node.find_farthest_node_with_overlap()
            else:
                node = next(iter(node.out_edges)).out_node
            if node is self.root:
                return
