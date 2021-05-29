""""""
from __future__ import annotations
from moPepGen import vep
from moPepGen.SeqFeature import FeatureLocation


class VariantRecordWithCoordinate():
    """"""
    def __init__(self, variant:vep.VEPVariantRecord, location:FeatureLocation):
        """"""
        self.variant = variant
        self.location = location

    def shift(self, index:int) -> VariantRecordWithCoordinate:
        """"""
        return self.__class__(
            variant=self.variant,
            location=FeatureLocation(
                start=self.location.start + index,
                end=self.location.end + index
            )
        )
    