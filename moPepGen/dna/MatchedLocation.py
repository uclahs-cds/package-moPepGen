""" MatchedLocation """
from __future__ import annotations
from moPepGen.SeqFeature import FeatureLocation


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
        start, stop, step = index.indices(len(self))
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