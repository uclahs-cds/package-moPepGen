""" Module for variant record with coordinate. """
from __future__ import annotations
from typing import TYPE_CHECKING
import math
from moPepGen.SeqFeature import FeatureLocation


if TYPE_CHECKING:
    from moPepGen.seqvar import VariantRecord

class VariantRecordWithCoordinate():
    """ This class models the variant record with its coordinate location at
    a gene or protein. This is used mainly in the graph to keep track on the
    location of variants of a node when the variable bubble expands forward
    or backward. """
    def __init__(self, variant:VariantRecord, location:FeatureLocation,
            is_stop_altering:bool=False, is_silent:bool=False,
            downstream_cleavage_altering:bool=False,
            upstream_cleavage_altering:bool=False):
        """ Constructor """
        self.variant = variant
        self.location = location
        self.is_stop_altering = is_stop_altering
        self.is_silent = is_silent
        self.downstream_cleavage_altering = downstream_cleavage_altering
        self.upstream_cleavage_altering = upstream_cleavage_altering

    def not_cleavage_altering(self) -> bool:
        """ Whether is not downstream or upstream cleavage altering """
        return not (self.upstream_cleavage_altering or self.downstream_cleavage_altering)

    def shift(self, index:int) -> VariantRecordWithCoordinate:
        """ Shift the coordinate of the object by a given number. """
        return self.__class__(
            variant=self.variant,
            location=FeatureLocation(
                start=max(self.location.start + index, 0),
                end=self.location.end + index,
                seqname=self.location.seqname,
                reading_frame_index=self.location.reading_frame_index,
                start_offset=self.location.start_offset,
                end_offset=self.location.end_offset
            ),
            is_stop_altering=self.is_stop_altering,
            is_silent=self.is_silent,
            upstream_cleavage_altering=self.upstream_cleavage_altering,
            downstream_cleavage_altering=self.downstream_cleavage_altering
        )

    def to_protein_coordinates(self) -> VariantRecordWithCoordinate:
        """ Returns a new object with the coordinates of the translated protein
        """
        start = math.floor(self.location.start / 3)
        end = math.ceil(self.location.end / 3)
        start_offset = self.location.start - start * 3
        end_offset = end * 3 - self.location.end
        return VariantRecordWithCoordinate(
            variant=self.variant,
            location=FeatureLocation(
                start=start, end=end, seqname=self.location.seqname,
                reading_frame_index=self.location.reading_frame_index,
                start_offset=start_offset, end_offset=end_offset
            ),
            is_stop_altering=self.is_stop_altering,
            is_silent=self.is_silent
        )

    def __getitem__(self, index) -> VariantRecordWithCoordinate:
        """ Get item """
        start, stop, _ = index.indices(self.location.end)
        start_index = max([start, self.location.start])
        end_index = min([stop, self.location.end])
        start_offset = self.location.start_offset if start == self.location.start else 0
        end_offset = self.location.end_offset if stop == self.location.end else 0
        location = FeatureLocation(
            start=start_index, end=end_index, seqname=self.location.seqname,
            reading_frame_index=self.location.reading_frame_index,
            start_offset=start_offset, end_offset=end_offset
        )
        return VariantRecordWithCoordinate(
            variant=self.variant,
            location=location,
            is_stop_altering=self.is_stop_altering,
            is_silent=self.is_silent
        )
