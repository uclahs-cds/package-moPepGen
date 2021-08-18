""" Module for SeqFeature """
from __future__ import annotations
from Bio.SeqFeature import FeatureLocation as BioFeatureLocation
from Bio.SeqFeature import SeqFeature as BioSeqFeature


_STRAND_LEVELS = {None:0, 0:1, -1:2, 1:3}

class FeatureLocation(BioFeatureLocation):
    """ This models the range of a sequence, with a start and an end. """
    def __init__(self, *args, seqname:str=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.seqname = seqname

    def __eq__(self, other: FeatureLocation)->bool:
        """ equal to """
        return self.start == other.start \
            and self.end == other.end \
            and self.strand == other.strand

    def __ne__(self, other: FeatureLocation) -> bool:
        """ Not equal to """
        return not self == other

    def __gt__(self, other: FeatureLocation)->bool:
        """ greater than """
        if self.start > other.start:
            return True
        if self.start == other.start:
            if _STRAND_LEVELS[self.strand] > _STRAND_LEVELS[other.strand]:
                return True
            if _STRAND_LEVELS[self.strand] == _STRAND_LEVELS[other.strand] \
                and self.end > other.end:
                return True
        return False

    def __ge__(self, other:FeatureLocation)->bool:
        """ greater or equal """
        return self > other or self == other

    def __lt__(self, other: FeatureLocation) -> bool:
        """ less than """
        return not self >= other

    def __le__(self, other: FeatureLocation) -> bool:
        """ less or equal to """
        return not self > other

    def __hash__(self):
        """ hash """
        return hash((self.seqname, self.start, self.end, self.strand))

    def overlaps(self, other:FeatureLocation) -> bool:
        """ Find whether the location overlaps with the other """
        return self.start in other or self.end in other or \
            other.start in self or other.end in self

    def get_overlap(self, other:FeatureLocation) -> FeatureLocation:
        """ Returns the range that is overlap """
        if not self.overlaps(other):
            raise ValueError('Two ranges do not overlap')
        return self.__class__(
            start=max(self.start, other.start),
            end=min(self.end, other.end)
        )

    def is_superset(self, other:FeatureLocation) -> bool:
        """ Find whether the location is a superset of the other """
        return self.start <= other.start and self.end >= other.end


class SeqFeature(BioSeqFeature):
    """ Models the annotation of a given range of the sequence. This extends
    the Bio.SeqFeature.SeqFeature object by adding a chromosome and some
    additional attributes. """
    def __init__(self, chrom:str, attributes:dict, *args, **kwargs):
        """ Constructor """
        self.location = None
        super().__init__(*args, **kwargs)
        if self.location is not None:
            self.location.__class__ = FeatureLocation
        self.location:FeatureLocation
        self.chrom = chrom
        self.attributes = attributes

    @property
    def biotype(self):
        """ biotype """
        return self.attributes['gene_type']

    def __eq__(self, other:SeqFeature) -> bool:
        """ equal to """
        return self.location == other.location

    def __ne__(self, other:SeqFeature) -> bool:
        """ not equal to """
        return not self.location == other.location

    def __gt__(self, other:SeqFeature) -> bool:
        """ greater than """
        return self.location > other.location

    def __ge__(self, other:SeqFeature) -> bool:
        """ greater or equal to """
        return self.location == other.location or \
            self.location > other.location

    def __lt__(self, other:SeqFeature) -> bool:
        """ less than """
        return not self >= other

    def __le__(self, other:SeqFeature) -> bool:
        """ less or equal to """
        return not self > other

    def __hash__(self):
        """ hash """
        return hash((self.chrom, self.location))
