""" CVGNode: node for ThreeFrameCVG """
from __future__ import annotations
import copy
from typing import List, Set
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen.seqvar import VariantRecord, VariantRecordWithCoordinate
from moPepGen.svgraph.TVGNode import TVGNode


class CVGNode(TVGNode):
    """ CVGNode """
    def __init__(self, circ:VariantRecord, seq:DNASeqRecordWithCoordinates,
            variants:List[VariantRecordWithCoordinate]=None,
            frameshifts:Set[VariantRecord]=None, branch:bool=False,
            orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None):
        """ constructor """
        super().__init__(
            seq=seq,
            variants=variants,
            frameshifts=frameshifts,
            orf=orf,
            reading_frame_index=reading_frame_index,
            subgraph_id=subgraph_id
        )
        self.circ = circ

    def is_reference(self) -> bool:
        """"""
        return not any(v.variant is not self.circ for v in self.variants)

    def copy(self) -> bool:
        """ copy """
        return self.__class__(
            circ=self.circ,
            seq=self.seq,
            variants=copy.copy(self.variants),
            frameshifts=copy.copy(self.frameshifts),
            branch=self.branch,
            orf=self.orf,
            reading_frame_index=self.reading_frame_index,
            subgraph_id=self.subgraph_id
        )

    def create_node(self, seq:DNASeqRecordWithCoordinates,
            variants:List[VariantRecordWithCoordinate]=None,
            frameshifts:Set[VariantRecord]=None, branch:bool=False,
            orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None):
        """ Constructor for TVGNode.

        Args:
            seq (DNASeqRecord): The sequence.
            variant (VariantRecord | None): The variant record or None for
                reference.
        """
        return self.__class__(
            circ=self.circ,
            seq=seq,
            variants=variants,
            frameshifts=frameshifts,
            branch=branch,
            orf=orf,
            reading_frame_index=reading_frame_index,
            subgraph_id=subgraph_id
        )

    def append_right(self, other:CVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq

        for variant in other.variants:
            if self.variants and \
                    self.variants[-1].location.end == variant.location.start:
                self.variants[-1].location = FeatureLocation(
                    start=self.variants[-1].location.start,
                    end=variant.location.end
                )
                continue
            self.variants.append(variant.shift(len(self.seq.seq)))

        self.seq = new_seq

    def append_left(self, other:TVGNode) -> None:
        """ Combine the other node the the left. """
        self.seq = other.seq + self.seq

        variants = copy.copy(other.variants)
        for variant in self.variants:
            if other.variants and \
                    variant.location.end == other.variants[-1].location.start:
                other.variants[-1].location = FeatureLocation(
                    start=other.variants[-1].location.start,
                    end=variant.location.end
                )
                continue
            variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants
