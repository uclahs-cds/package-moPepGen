""""""
from __future__ import annotations
from moPepGen import seqvar
from moPepGen.SeqFeature import FeatureLocation


class VariantRecordWithCoordinate():
    """"""
    def __init__(self, variant:seqvar.VariantRecord, location:FeatureLocation):
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
    