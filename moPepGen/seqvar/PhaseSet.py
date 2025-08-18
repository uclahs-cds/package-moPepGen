""" Module for handling phase sets in variant records """
from __future__ import annotations
from typing import TYPE_CHECKING
from moPepGen import err


if TYPE_CHECKING:
    from typing import Set, Iterable

def have_competible_phase_sets(set1:Set[str], set2:Set[str],
        groups:Iterable[Set[str]]) -> bool:
    """ Check if two sets of variants have compatible phase sets. """
    for group in groups:
        s1 = set1 & group
        s2 = set2 & group
        if s1 and s2 and not (s1 & s2):
            return False
    return True

def join_phase_sets(set1:Set[str], set2:Set[str], group:Iterable[Set[str]]) -> Set[str]:
    """ Join two sets of phase sets """
    joined = set()
    for group in groups:
        s1 = set1 & group
        s2 = set2 & group
        if s1 and not s2:
            joined.update(s1)
        elif s2 and not s1:
            joined.update(s2)
        elif s1 and s2:
            s = s1 & s2
            if not s:
                raise err.PhaseSetsUncompatibleError(s1, s2)
            joined.update(s)
    return joined
