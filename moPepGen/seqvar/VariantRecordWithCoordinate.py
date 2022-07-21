""" Module for variant record with coordinate. """
from __future__ import annotations
import math
from moPepGen import seqvar
from moPepGen.SeqFeature import FeatureLocation


class VariantRecordWithCoordinate():
    """ This class models the variant record and it's coordinate location at
    a sequence. This is used mainly in the graph, to keep track on the location
    of variants of a node when it expand forward or backward. """
    def __init__(self, variant:seqvar.VariantRecord, location:FeatureLocation,
            is_stop_altering:bool=False):
        """ Constructor """
        self.variant = variant
        self.location = location
        self.is_stop_altering = is_stop_altering

    def shift(self, index:int) -> VariantRecordWithCoordinate:
        """ Shift the coordinate of the object by a given number. """
        return self.__class__(
            variant=self.variant,
            location=FeatureLocation(
                start=max(self.location.start + index, 0),
                end=self.location.end + index
            ),
            is_stop_altering=self.is_stop_altering
        )

    def to_protein_coordinates(self) -> VariantRecordWithCoordinate:
        """ Returns a new object with the coordinates of the translated protein
        """
        start = int(self.location.start / 3)
        end = math.ceil(self.location.end / 3)
        return seqvar.VariantRecordWithCoordinate(
            variant=self.variant,
            location=FeatureLocation(start=start, end=end),
            is_stop_altering=self.is_stop_altering
        )

    def __getitem__(self, index) -> VariantRecordWithCoordinate:
        """ Get item """
        start, stop, _ = index.indices(self.location.end)
        start_index = max([start, self.location.start])
        end_index = min([stop, self.location.end])
        location = FeatureLocation(start=start_index, end=end_index)
        return VariantRecordWithCoordinate(
            variant=self.variant,
            location=location
        )
