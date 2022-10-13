""" """
from __future__ import annotations
from typing import Set, List
import copy
from moPepGen import seqvar


class PVGOrf():
    """ """
    def __init__(self, orf:List[int,int]=None,
            start_gain:Set[seqvar.VariantRecord]=None):
        """ constructor """
        self.orf = orf or None
        self.start_gain = start_gain or set()

    def copy(self) -> PVGOrf:
        """ copy """
        orf = self.orf
        start_gain = copy.copy(self.start_gain)
        return self.__class__(orf, start_gain)

    def __eq__(self, other:PVGOrf) -> bool:
        """ equal """
        return self.orf == other.orf and self.start_gain == other.start_gain

    def __ne__(self, other:PVGOrf) -> bool:
        """ not equal """
        return not self == other

    def __gt__(self, other:PVGOrf) -> bool:
        """ Greater than. We want it to be "greater" if the ORF is larger and
        less start gain variants or larger variant position. """
        if self.orf[0] > other.orf[0]:
            return True
        if self.orf[0] < other.orf[0]:
            return False
        if len(self.start_gain) < len(other.start_gain):
            return True
        if len(self.start_gain) > len(other.start_gain):
            return False
        for x, y in zip(self.start_gain, other.start_gain):
            if x > y:
                return True
            if x < y:
                return False
        return False

    def __ge__(self, other:PVGOrf) -> bool:
        """ larger then or equal to """
        return self > other or self == other

    def __lt__(self, other:PVGOrf) -> bool:
        """ less than """
        return not self >= other

    def __le__(self, other:PVGOrf) -> bool:
        """ less then or equal to """
        return not self > other

    def __hash__(self):
        """ hash """
        return hash((tuple(self.orf), frozenset(self.start_gain)))
