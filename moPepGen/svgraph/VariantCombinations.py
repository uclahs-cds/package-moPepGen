""" Frameshifting mutations of a node """
from __future__ import annotations
import copy
from typing import FrozenSet, List, Set
from moPepGen import seqvar


class VariantCombinations():
    """ Node Framesfhits """
    def __init__(self, data:Set[FrozenSet[seqvar.VariantRecord]]=None, max_n:int=3):
        """ """
        self.data = data or set()
        self.max_n = max_n

    def add(self, variant:seqvar.VariantRecord) -> None:
        """ Combine a variant to each of the variant combinations """
        data = set()
        for comb in self.data:
            if len(comb) < self.max_n:
                comb = set(comb)
                comb.add(variant)
                comb = frozenset(comb)
            data.add(comb)
        self.data = data

    def extend(self, variants:VariantCombinations) -> None:
        """ Extend another VariantCombinations  """
        new_one = copy.copy(variants)
        self.data.update(new_one.data)

    def __copy__(self) -> VariantCombinations:
        """ create a copy """
        data = set()
        for comb in self.data:
            data.add(copy.copy(comb))
        return self.__class__(data)

    def remove(self, variant:seqvar.VariantRecord):
        """ remove """
        data = set()
        for comb in self.data:
            if variant in comb:
                comb = set(comb)
                comb.remove(variant)
                if len(comb) > 0:
                    comb = frozenset(comb)
                    data.add(comb)
            else:
                data.add(comb)
        self.data = data

    def update(self, variants:List[seqvar.VariantRecordWithCoordinate],
            add_singleton:bool=True) -> None:
        """ For a given VariantRecordWithCoordinates, update those that are
        frameshift mutations """
        for variant in variants:
            if variant.variant.is_frameshifting():
                self.add(variant.variant)
            if add_singleton:
                singleton = frozenset([variant.variant])
                if singleton not in self.data:
                    self.data.add(singleton)
