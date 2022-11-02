""" Module for CircularVariantGraph, a directed cyclic graph, for circRNA etc.
"""
from __future__ import annotations
import copy
from typing import Dict, Union, List, TYPE_CHECKING
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import svgraph, seqvar, params
from moPepGen.svgraph.TVGNode import TVGNode
from moPepGen.svgraph.SubgraphTree import SubgraphTree


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
            circ_record:circ.CircRNAModel=None, attrs:dict=None,
            subgraphs:SubgraphTree=None, cleavage_params:params.CleavageParams=None):
        """ Construct a CircularVariantGraph

        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the
                    transcript (reference).
            transcript_id (str): The transcript ID that the circRNA is
                associated with.
            circ (CircRNAodel): The circRNA model.
            attrs (dict): additional attributes

        """
        self.seq = seq
        self.id = _id
        if self.seq and not self.seq.locations:
            self.add_default_sequence_locations()
        self.attrs = attrs
        self.circ = circ_record
        location = FeatureLocation(
            seqname=self.id,
            start=seq.locations[0].ref.start,
            end=seq.locations[0].ref.end
        )
        circ_variant = seqvar.VariantRecord(
            location=location,
            ref=seq.seq[0],
            alt='<circRNA>',
            _type='circRNA',
            _id=_id
        )
        super().__init__(
            seq=seq, _id=_id, root=root, reading_frames=reading_frames,
            cds_start_nf=cds_start_nf, has_known_orf=has_known_orf,
            global_variant=circ_variant, subgraphs=subgraphs,
            cleavage_params=cleavage_params
        )

    def get_circ_variant_with_coordinate(self) -> seqvar.VariantRecordWithCoordinate:
        """ Add a variant record to the frameshifting of the root node. This
        will treat all peptides as variant peptide in the later steps. """
        location = FeatureLocation(
            seqname=self.circ.gene_id,
            start=0,
            end=len(self.seq)
        )
        return seqvar.VariantRecordWithCoordinate(self.global_variant, location)

    def init_three_frames(self, truncate_head:bool=False):
        """ Initiate the three reading-frame graph. """
        super().init_three_frames(truncate_head=truncate_head)
        var = self.get_circ_variant_with_coordinate()
        for root in self.reading_frames:
            if len(root.out_edges) > 1:
                raise ValueError('Initiating CVG should not contain any variant')
            node = list(root.out_edges)[0].out_node
            node.variants.append(var)
            self.add_edge(node, root, 'reference')

    def create_variant_circ_graph(self, variants:List[seqvar.VariantRecord]):
        """ Apply a list of variants to the graph. Variants not in the
        range are ignored. Variants at the first nucleotide of each fragment
        of the sequence are also ignored, because it causes the exon splice
        site to be changed.
        """
        filtered_variants = []
        exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion', 'circRNA']
        for variant in variants:
            if variant.type == exclude_type:
                continue

            circ_start = self.circ.fragments[0].location.start
            if circ_start <= variant.location.start < circ_start + 3:
                continue

            for fragment in self.circ.fragments:
                location = fragment.location
                if location.start <= variant.location.start < location.end:
                    filtered_variants.append(variant)
                    break
        super().create_variant_graph(filtered_variants, None, None, None)

    def get_tail_nodes(self) -> Dict[TVGNode, TVGNode]:
        """ The the last node of each reading frame. """
        tail_nodes:List[TVGNode] = [None, None, None]
        for i, root in enumerate(self.reading_frames):
            if len(root.in_edges) > 2:
                raise ValueError('CVG should not have any stop altering mutation')
            in_edge = None
            for in_edge in root.in_edges:
                if in_edge.in_node is not self.root:
                    break
            if not in_edge:
                raise ValueError('The TVG is not cyclic')
            tail_nodes[i] = in_edge.in_node
            self.remove_edge(in_edge)
        return tail_nodes

    def extend_loop(self):
        """ Extend the three frame circular graph into linear by repeating the
        graph of each reading frame for 3 additional copies. """
        n_loops = 3

        subgraphs:List[ThreeFrameCVG] = []
        cur = self
        for _ in range(n_loops):
            sub = cur.create_subgraph_copy()
            subgraphs.append(sub)
            cur = sub

        cur = self
        cur_tail_nodes = self.get_tail_nodes()
        cur = self
        for sub in subgraphs:
            sub_tail_nodes = sub.get_tail_nodes()
            self.subgraphs.add_subgraph(
                child_id=sub.id, parent_id=cur.id, level=sub.root.level,
                start=self.seq.locations[-1].ref.end, end=self.seq.locations[-1].ref.end
            )
            for edge in sub.root.out_edges:
                root = edge.out_node
                if len(root.out_edges) > 1:
                    raise ValueError('CVG should not have any start altering mutation')
                head = list(root.out_edges)[0].out_node
                end_node = cur_tail_nodes[head.reading_frame_index]
                for _edge in copy.copy(head.in_edges):
                    cur.remove_edge(_edge)
                cur.add_edge(end_node, head, 'reference')
            cur = sub
            cur_tail_nodes = sub_tail_nodes

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

    def create_subgraph_copy(self) -> ThreeFrameCVG:
        """ Create a subgraph by copying the current graph. The subgraph created
        will be attached to the end of this current graph. """
        subgraph_id = self.subgraphs.generate_subgraph_id()
        seq = copy.deepcopy(self.seq)
        for loc in seq.locations:
            loc.ref.seqname = subgraph_id
        root = self.root.deepcopy(subgraph_id, 1)
        reading_frames:Dict[TVGNode] = [None, None, None]
        for node in root.get_out_nodes():
            reading_frames[node.reading_frame_index] = node
        return ThreeFrameCVG(
            seq=seq, _id=subgraph_id, root=root, reading_frames=reading_frames,
            cds_start_nf=self.cds_start_nf, has_known_orf=self.has_known_orf,
            circ_record=self.circ, attrs=copy.copy(self.attrs),
            subgraphs=None, cleavage_params=self.cleavage_params
        )
