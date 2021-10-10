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

    def __str__(self) -> str:
        """ str """
        return f"{self.seqname}:{int(self.start)}-{int(self.end)}"

    def __hash__(self):
        """ hash """
        return hash((self.seqname, self.start, self.end, self.strand))

    def overlaps(self, other:FeatureLocation) -> bool:
        """ Find whether the location overlaps with the other """
        return self.start in other or self.end - 1 in other or \
            other.start in self or other.end - 1 in self

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

    def _shift(self, offset:int) -> SeqFeature:
        """ shift by i """
        new_feature = super()._shift(offset)
        new_feature.__class__ = self.__class__
        return new_feature

class MatchedLocation():
    """ A MappedLocation object holds the location of a sequence that is
    matched with another. It is mainly used in DNASeqRecordWithCoordinates
    class to represent the part of the sequence that matches with the
    reference sequence.

    For example, a sequence of length 20, mapped to position 1917 to 1937 of
    the transcript, can be represented as:

    query = FeatureLocation(start=0, end=20)
    ref = FeatureLocation(start=1917, end=1937)
    location = MatchedLocation(query=query, ref=ref)

    Attributes:
        query (FeatureLocation): The location of the query sequence.
        ref (FeatureLocation): The location of the reference sequence.
    """
    def __init__(self, query:FeatureLocation, ref:FeatureLocation):
        """ Constructor for MatchedLocation

        Args:
            query (FeatureLocation): The location of the query sequence.
            ref (FeatureLocation): The location of the reference sequence.
        """
        if len(query) != len(ref):
            raise ValueError('Location length must equal.')
        self.query = query
        self.ref = ref

    def __len__(self):
        """ length """
        return len(self.query)

    def __eq__(self, other:MatchedLocation) -> bool:
        """ equal to """
        return self.ref == other.ref

    def __gt__(self, other:MatchedLocation) -> bool:
        """ greater than """
        return self.ref > other.ref

    def __getitem__(self, index) -> MatchedLocation:
        """ Get item. The query location of the returned object starts at 0.
        """
        start, stop, _ = index.indices(len(self))
        query = FeatureLocation(
            seqname=self.query.seqname,
            start=0,
            end=stop - start
        )
        ref = FeatureLocation(
            seqname=self.ref.seqname,
            start=self.ref.start + start,
            end=self.ref.start + stop
        )
        return self.__class__(
            query=query,
            ref=ref
        )

    def shift(self, i:int) -> MatchedLocation:
        """ Shift query window by i """
        query = self.query.__class__(
            seqname=self.query.seqname,
            start=self.query.start._shift(i),
            end=self.query.end._shift(i)
        )
        return self.__class__(query=query, ref=self.ref)
