""" Module for CircularVariantGraph, a directed cyclic graph, for circRNA etc.
"""
from __future__ import annotations
import copy
from typing import Dict, Union, List, TYPE_CHECKING
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import svgraph, seqvar
from moPepGen.svgraph.TVGNode import TVGNode


if TYPE_CHECKING:
    from moPepGen import circ
    from moPepGen.dna import DNASeqRecordWithCoordinates

class ThreeFrameCVG(svgraph.ThreeFrameTVG):
    """ Defines a directed cyclic graph for circular nucleotide molecules such
    ass circRNA and the variants associated with it.

    Attributes:
        seq (DNASeqRecordWithCoordinates): The original sequence of the
                transcript (reference).
        transcript_id (str): The transcript ID that the circRNA is associated
            with
        attrs (dict): additional attributes
    """
    def __init__(self, seq:Union[DNASeqRecordWithCoordinates,None],
            _id:str, root:TVGNode=None, reading_frames:List[TVGNode]=None,
            cds_start_nf:bool=False, has_known_orf:bool=False,
            circ_record:circ.CircRNAModel=None, attrs:dict=None):
        """ Construct a CircularVariantGraph

        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the
                    transcript (reference).
            transcript_id (str): The transcript ID that the circRNA is
                associated with.
            circ (CircRNAodel): The circRNA model.
            attrs (dict): additional attributes

        """
        super().__init__(seq, _id, root, reading_frames, cds_start_nf, has_known_orf)
        self.attrs = attrs
        self.circ = circ_record

    def init_three_frames(self, truncate_head:bool=False):
        super().init_three_frames(truncate_head)
        var = self.as_frameshifting()
        for root in self.reading_frames:
            if len(root.out_edges) > 1:
                raise ValueError('Initiating CVG should not contain any variant')
            node = list(root.out_edges)[0].out_node
            node.variants.append(var)
            self.add_edge(node, root, 'reference')

    def as_frameshifting(self) -> seqvar.VariantRecordWithCoordinate:
        """ Add a variant record to the frameshifting of the root node. This
        will treat all peptides as variant peptide in the later steps. """
        location = FeatureLocation(
            seqname=self.circ.gene_id,
            start=self.seq.locations[0].ref.start,
            end=self.seq.locations[0].ref.end
        )
        variant = seqvar.VariantRecord(
            location=location,
            ref=self.seq.seq[0],
            alt='<circRNA>',
            _type='circRNA',
            _id=self.id
        )
        location = FeatureLocation(
            seqname=self.circ.gene_id,
            start=0,
            end=len(self.seq)
        )
        return seqvar.VariantRecordWithCoordinate(variant, location)

    def create_variant_circ_graph(self, variants: List[seqvar.VariantRecord]):
        """ Apply a list of variants to the graph. Variants not in the
        range are ignored. Variants at the first nucleotide of each fragment
        of the sequence are also ignored, because it causes the exon splice
        site to be changed.
        """
        filtered_variants = []
        for variant in variants:
            if variant.type == 'Fusion':
                continue
            for location in self.seq.locations:
                if variant.location.start > location.ref.start + 3 and \
                        variant.location.end < location.ref.end:
                    filtered_variants.append(variant)
                    break
        super().create_variant_graph(filtered_variants, None, None, None)

    def extend_loop(self):
        """ Extend each reading frame for one more loop. """
        frame_map:Dict[TVGNode, TVGNode] = {}
        for root in self.reading_frames:
            if len(root.in_edges) > 2:
                raise ValueError('CVG should not have any stop altering mutation')
            in_edge = None
            for in_edge in root.in_edges:
                if in_edge.in_node is not self.root:
                    break
            if not in_edge:
                raise ValueError('The TVG is not cyclic')
            frame_map[root] = in_edge.in_node
            self.remove_edge(in_edge)

        root_copy = self.root.deepcopy()
        for edge in root_copy.out_edges:
            root = edge.out_node
            if len(root.out_edges) > 1:
                raise ValueError('CVG should not have any start altering mutation')
            head = list(root.out_edges)[0].out_node
            end_node = frame_map[self.reading_frames[head.reading_frame_index]]
            for _edge in copy.copy(head.in_edges):
                self.remove_edge(_edge)
            self.add_edge(end_node, head, 'reference')

    def truncate_three_frames(self):
        """ For each of the reading frames, tuncate the leading nodes for 1
        nucleotide. """
        for i, root in enumerate(self.reading_frames):
            if i == 0:
                continue
            if len(root.out_edges) > 1:
                raise ValueError('CVG should not contain any start altering mutaiton')
            head = list(root.out_edges)[0].out_node
            head.truncate_left(i)
