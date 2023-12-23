""" DNASeqRecord """
from __future__ import annotations
import copy
from typing import Iterable, List
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen.aa.expasy_rules import EXPASY_RULES
from moPepGen import aa


class DNASeqRecord(SeqRecord):
    #pylint: disable=W0223
    """ A DNASeqRecord object holds a single DNA sequence and information about
    it. Inherits the Bio.SeqRecord.SeqRecord class.

    Attributes:
        seq (Bio.Seq.Seq: The DNA sequence.
        id (str)
        name (str)
        description (str)

    For other attributes see biopython's documentation:
    https://biopython.org/docs/1.75/api/Bio.Seq.html?highlight=seq#Bio.Seq.Seq
    """

    def __add__(self, other:DNASeqRecord)->DNASeqRecord:
        """ add operation """
        new = super().__add__(other)
        new.__class__ = DNASeqRecord
        return new

    def __hash__(self):
        """ hash """
        return hash(str(self.seq))

    def reverse_complement(self, _id=False, name=False, description=False,
            features=True, annotations=False, letter_annotations=True,
            dbxrefs=False)->DNASeqRecord:
        """ Returns the reverse complement sequence """
        seq = super().reverse_complement(
            id=_id, name=name, description=description, features=features,
            annotations=annotations, letter_annotations=letter_annotations,
            dbxrefs=dbxrefs
        )
        seq.__class__ = self.__class__
        return seq

    def find_start_codon(self) -> int:
        """ find start codon index """
        return self.seq.find('ATG')

    def iter_start_codon(self) -> Iterable[int]:
        """ Create iterators of start codon positions. """
        for it in re.finditer('ATG', str(self.seq)):
            yield it.start()

    def find_stop_codon(self, start:int=0) -> int:
        """ Find and return the stop codon.

        Args:
            start (int): The start position to search. It only searches the
                reading frame specified by start.
            table (str): The codon table to choose (unsupported)

        Returns:
            The index of the stop codon. -1 is returned if not found.
        """
        # NOTE(CZ): enable choosing codon table from the Bio.Data.CodonTable
        # module.
        stop_codons = ['TAA', 'TAG', 'TGA']
        n = len(self)
        for i in range(start, n - n % 3, 3):
            if str(self.seq[i:i+3]) in stop_codons:
                return i
        return -1

    def find_orf_first(self, start:int=0) -> int:
        """ Search for the first open reading frame (start codon)

        Args:
            start (int): The starting index to search.

        Returns:
            The index of the start codon. -1 is returned if not found.
        """
        # NOTE(CZ): consider other start codons?
        seq:Seq = self.seq
        return seq.find('ATG', start=start)

    def find_orf_all(self, start:int=0) -> int:
        """ Returns iterator of start codon positions """

        while True:
            start = self.find_orf_first(start)
            if start == -1:
                return
            yield start
            start += 3

    def iter_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            start:int=None, end:int=None) -> Iterable[int]:
        """ Returns a generator of the enzymatic cleave sites of a given range.
        Default is the entire sequence.

        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            start (int): Index to start searching.
            end (int): Index to stop searching.

        Returns:
            A generator of cleave sites.
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self)
        peptides = str(self.seq[start:end].translate())
        rule = EXPASY_RULES[rule]
        exception = EXPASY_RULES.get(exception, exception)
        exceptions = [] if exception is None else \
            [x.end() for x in re.finditer(exception, peptides)]

        for x in re.finditer(rule, peptides):
            if x not in exceptions:
                yield x.end() * 3 + start

    def find_all_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            start:int=None, end:int=None) -> List[int]:
        """ Find all enzymatic cleave sites. """
        if start is None:
            start = 0
        if end is None:
            end = len(self)
        end = end - (end - start) % 3
        return list(self.iter_enzymatic_cleave_sites(rule=rule,
            exception=exception, start=start, end=end))

    def find_all_start_codons(self) -> List[int]:
        """ Find all start codon positions """
        start_positions = []
        start_codon = 'ATG'
        while True:
            if len(start_positions) == 0:
                i = 0
            else:
                i = start_positions[-1] + 1
            j = self.seq[i:].find(start_codon)
            if j == -1:
                break
            start_positions.append(i + j)
        return start_positions

    def find_last_cleave_position(self, end:int, rule:str, exception:str=None,
            miscleavage:int=0) -> int:
        """ Find the last enzymatic cleave site with a given number of
        miscleavage. """
        end = end - end % 3
        cleave_sites = [0]
        cleave_sites.extend(self.find_all_enzymatic_cleave_sites(rule=rule,
            exception=exception, start=None, end=end))
        i = len(cleave_sites) - miscleavage
        return max(i, 0)

    def find_next_cleave_position(self, start:int, rule:str,
            exception:str=None, miscleavage:int=0) -> int:
        """ Find the next enzymatic cleave site with a given number of
        miscleavage. Start must be in the correct ORF to return the correct
        site. """
        i = 0
        for site in self.iter_enzymatic_cleave_sites(rule=rule,
                exception=exception, start=start, end=None):
            if i == miscleavage:
                return site + start
            i += 1
        return -1

    def find_cleave_positions_before(self, position:int, rule:str,
            exception:str=None, miscleavage:str=2) -> List[int]:
        """ Find specified number of enzymatic cleave sites before a given
        position.
        """
        sites = []
        position = position - position % 3
        left_cleave_sites = [0]
        left_cleave_sites.extend(self.find_all_enzymatic_cleave_sites(
            rule=rule, exception=exception, start=None, end=position))
        i = miscleavage +1
        while i > 0:
            if len(left_cleave_sites) >= i:
                sites.append(left_cleave_sites[-i])
            else:
                sites.append(None)
            i += 1
        return sites

    def find_cleave_positions_after(self, position:int, rule:str,
            exception:str=None, miscleavage:str=2) -> List[int]:
        """ Find specified number of enzymatic cleave sites before a given
        position. """
        sites = []
        start = position - position % 3
        right_cleave_sites = self.find_all_enzymatic_cleave_sites(rule=rule,
            exception=exception, start=start, end=None)
        for i, _ in enumerate(sites):
            sites[i] += start
        right_cleave_sites.append(len(self))
        for i in range(miscleavage):
            if len(right_cleave_sites) >= i + 1:
                sites.append(right_cleave_sites[i])
            else:
                sites.append(None)
        return sites


    def find_cleave_positions_within(self, start:int, end:int, rule:int,
            exception:str=None) -> List[int]:
        """ Find enzymatic cleave sites within a given range """
        start = start - start % 3
        end = end + 3 - end % 3
        left = self.find_last_cleave_position(end=start, rule=rule,
            exception=exception, miscleavage=0)
        right = self.find_next_cleave_position(start=end, rule=rule,
            exception=exception, miscleavage=0)
        new_sites = self[left:right].find_all_enzymatic_cleave_sites(rule=rule,
            exception=exception)
        for i, _ in enumerate(new_sites):
            new_sites[i] += start
        return new_sites

    def translate(self, table:str='Standard', stop_symbol:str="*",
            to_stop:bool=False, cds:bool=False, gap:str="-"
            ) -> aa.AminoAcidSeqRecord:
        """ Translate the DNA sequence to a amino acid sequence.

        Args:
            table (str): Which codon table to use? This can be either a name
                (string), an NCBI identifier (integer), or a CodonTable object
                (useful for non-standard genetic codes). This defaults to the
                "Standard" table.
            to_stop (bool):  Boolean, defaults to False meaning do a full
                translation continuing on past any stop codons (translated as
                the specified stop_symbol). If True, translation is terminated
                at the first in frame stop codon (and the stop_symbol is not
                appended to the returned protein sequence).
        """
        if len(self.seq) % 3 > 0:
            end = len(self.seq) - len(self.seq) % 3
            seq:Seq = self.seq[:end]
            peptides = seq.translate(
                table=table, stop_symbol=stop_symbol, to_stop=to_stop,
                cds=cds, gap=gap
            )
        else:
            peptides = self.seq.translate(
                table=table, stop_symbol=stop_symbol, to_stop=to_stop,
                cds=cds, gap=gap
            )
        return aa.AminoAcidSeqRecord(
            seq=peptides,
            _id=self.id,
            name=self.name,
            description=self.description
        )


class DNASeqRecordWithCoordinates(DNASeqRecord):
    """ The DNASeqReocrdWithLocation holds a sequence and its location matched
    to another sequence, e.g., chromosome or transcript. It contains two
    additional attributes, locations and orf.

    Attributes:
        seq (Seq): The DNA sequence.
        locations (List[MatchedLocation]): The locations where this
            sequence is aligned tp.
        orf (FeatureLocation): The open reading frame start and end.
    """
    def __init__(self, seq:Seq, *args, locations:List[MatchedLocation]=None,
            orf:FeatureLocation=None, selenocysteine:List[FeatureLocation]=None,
            **kwargs,):
        """ Constract a DNASeqRecordWithCoordinates object.

        Args:
            seq (Seq): The DNA sequence.
            locations (List[MatchedLocation]): The locations where this
                sequence is aligned tp.
            orf (FeatureLocation): The open reading frame start and end.
        """
        super().__init__(seq=seq, *args, **kwargs)
        self.locations = locations or []
        # query index
        self.orf = orf
        self.selenocysteine = selenocysteine or []

    def __add__(self, other:DNASeqRecordWithCoordinates
            ) -> DNASeqRecordWithCoordinates:
        """ add operation """
        new = super().__add__(other)
        new.__class__ = DNASeqRecordWithCoordinates
        left_locs = copy.copy(self.locations)
        right_locs = copy.copy(other.locations)
        if left_locs and right_locs:
            lhs = self.locations[-1]
            rhs = other.locations[0].shift(len(self))
            if lhs.ref.end == rhs.ref.start and lhs.query.end == rhs.query.start \
                    and lhs.ref.seqname == rhs.ref.seqname \
                    and lhs.query.reading_frame_index == rhs.query.reading_frame_index:
                query = FeatureLocation(
                    start=lhs.query.start, end=rhs.query.end,
                    seqname=lhs.query.seqname,
                    reading_frame_index=lhs.query.reading_frame_index
                )
                ref = FeatureLocation(start=lhs.ref.start, end=rhs.ref.end,
                    seqname=lhs.ref.seqname)
                new_loc = MatchedLocation(query=query, ref=ref)
                right_locs.pop(0)
                left_locs[-1] = new_loc

        new.locations = left_locs + [loc.shift(len(self)) for loc in right_locs]
        new.orf = self.orf
        new.selenocysteine = self.selenocysteine
        return new

    def __getitem__(self, index)->DNASeqRecordWithCoordinates:
        """ get item """
        if isinstance(index, int):
            return super().__getitem__(index)
        start, stop, _ = index.indices(len(self))

        locations = []
        if start != stop:
            for location in self.locations:
                lhs = location.query.start
                rhs = location.query.end
                if rhs <= start:
                    continue
                if lhs >= stop:
                    break
                if lhs <= start:
                    if rhs <= stop:
                        location = location[start-lhs:]
                    else:
                        location = location[start-lhs:stop-lhs]
                else:
                    if rhs <= stop:
                        location = location.shift(-start)
                    else:
                        new_location = location[:stop-lhs]
                        new_location = new_location.shift(lhs-start)
                        location = new_location

                locations.append(location)

        return self.__class__(
            seq=self.seq[index],
            locations=locations,
            orf = self.orf,
            selenocysteine=self.selenocysteine,
            id=self.id,
            name=self.name,
            description=self.description
        )

    def __hash__(self):
        """hash"""
        return hash(self.seq)

    def __eq__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ equal to """
        if self.seq != other.seq:
            return False
        if len(self.locations) != len(other.locations):
            return False
        for i, j in zip(self.locations, other.locations):
            if i != j:
                return False
        return True

    def __ne__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ not equal to """
        return not self == other

    def __gt__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ greater than """
        for i, j in zip(self.locations, other.locations):
            if i > j:
                return True
            if i < j:
                return False
        if len(self.locations) > len(other.locations):
            return True
        return False

    def __ge__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ greater or equal to """
        return self > other or self == other

    def __lt__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ less than """
        return not self >= other

    def __le__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ less or equal """
        return not self > other

    def __le___(self, other:DNASeqRecordWithCoordinates):
        """ Due to an typo in biopython """
        return not self >= other

    def concat(self, other:DNASeqRecordWithCoordinates, _id:str=None,
            name:str=None, description:str=None
            ) -> DNASeqRecordWithCoordinates:
        """ concat two sequences into one """
        return self.__class__(
            seq=self.seq + other.seq,
            locations=self.locations + other.locations,
            orf=self.orf,
            id=self.id if _id is None else _id,
            name=self.name if name is None else name,
            description=self.description if description is None else \
                description
        )

    def search_orf(self):
        """ Searches the the open reading frame, finds the start and end
        position at the reference sequence and replace the `self.orf`. If
        no orf is found, (-1, -1) is returned. """
        start = self.find_orf_first()
        if start == -1:
            end = -1
        else:
            end = self.find_stop_codon(start)
            if end == -1:
                end = len(self) - len(self) % 3
        self.orf = FeatureLocation(
            start=start,
            end=end
        )

    def get_query_index(self, ref_index:int) -> int:
        """ Returns the query index wiht a given reference index """
        for location in self.locations:
            if ref_index in location.ref:
                return location.query.start + ref_index - location.ref.start
        return -1
