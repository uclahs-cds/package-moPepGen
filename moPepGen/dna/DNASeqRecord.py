""" DNASeqRecord """
from __future__ import annotations
from typing import Iterable, List, Tuple
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.aa.expasy_rules import EXPASY_RULES
from moPepGen import aa, dna


class DNASeqRecord(SeqRecord):
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
    
    def reverse_complement(self)->DNASeqRecord:
        """ Returns the reverse complement sequence """
        """ Get reverse complement sequence """
        seq = super().reverse_complement()
        seq.__class__ = self.__class__
        return seq
    
    def find_stop_codon(self, start:int=0, table:str='Standard') -> int:
        """ Find and return teh stop codon.
        
        Args:
            start (int): The start position to search. It only searches the
                reading frame specified by start.
            table (str): The codon table to choose (unsupported)

        Returns:
            The index of the stop codon. -1 is returned if not found.
        """
        # TODO(CZ): enable choosing codon table from the Bio.Data.CodonTable
        # module.
        stop_codons = ['TAA', 'TAG', 'TGA']
        n = len(self)
        for i in range(start, n - n % 3, 3):
            if self.seq[i:i+3] in stop_codons:
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
        return self.seq.find('ATG', start=start)
    
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
        return
    
    def find_all_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            start:int=None, end:int=None) -> List[int]:
        """ Find all enzymatic cleave sites. """
        if start is None:
            start = 0
        if end is None:
            end = len(self)
        end = end - (end - start) % 3
        return [i for i in self.iter_enzymatic_cleave_sites(rule=rule,
            exception=exception, start=start, end=end)]
    
    def find_last_cleave_position(self, end:int, rule:str, exception:str=None,
            miscleavage:int=0) -> int:
        """ Find the last enzymatic cleave site with a given number of
        miscleavage. """
        end = end - end % 3
        cleave_sites = [0]
        cleave_sites.extend(self.find_all_enzymatic_cleave_sites(rule=rule,
            exception=exception, start=None, end=end))

        i = len(cleave_sites) - miscleavage
        if i < 0:
            i = 0
        return i
    
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
        for i in range(len(sites)):
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
        for i in range(len(new_sites)):
            new_sites[i] += start
        return new_sites
    
    def translate(self, table:str='Standard', to_stop:bool=False, id:str=None,
            name:str=None, description:str=None, protein_id:str=None,
            transcript_id:str=None, gene_id:str=None) -> aa.AminoAcidSeqRecord:
        """"""
        if len(self.seq) % 3 > 0:
            end = len(self.seq) - len(self.seq) % 3
            peptides = self.seq[:end].translate(table=table, to_stop=to_stop)
        else:
            peptides = self.seq.translate(table=table, to_stop=to_stop)
        return aa.AminoAcidSeqRecord(
            seq=peptides,
            id=id if id else self.id,
            name=name if name else self.name,
            description=description if description else self.description,
            protein_id=protein_id,
            transcript_id=transcript_id,
            gene_id=gene_id
        )


class DNASeqRecordWithCoordinates(DNASeqRecord):
    """ The DNASeqReocrdWithLocation holds a sequence and its location matched
    to another sequence, e.g., chromosome or transcript.
    """
    def __init__(self, seq:Seq, locations:List[dna.MatchedLocation],
        orf:FeatureLocation=None, *args, **kwargs):
        """  """
        super().__init__(seq=seq, *args, **kwargs)
        self.locations = locations
        # query index
        self.orf = orf

    def __add__(self, other:DNASeqRecordWithCoordinates
            ) -> DNASeqRecordWithCoordinates:
        """ add operation """
        new = super().__add__(other)
        new.__class__ = DNASeqRecordWithCoordinates
        new.locations = self.locations + \
            [location.shift(len(self)) for location in other.locations]
        new.orf = self.orf
        return new
    
    def __getitem__(self, index)->DNASeqRecordWithCoordinates:
        """ get item """
        if isinstance(index, int):
            return super().__getitem__(index)
        start, stop, step = index.indices(len(self))

        locations = []
        if start != stop:
            for location in self.locations:
                lhs = location.query.start
                rhs = location.query.end
                if rhs < start:
                    continue
                if lhs > stop:
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

    def __gt__(self, other:DNASeqRecordWithCoordinates) -> bool:
        """ greater than """
        for i, j in zip(self.locations, other.locations):
            if i > j:
                return True
            elif i < j:
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
    
    def concat(self, other:DNASeqRecordWithCoordinates, id:str=None,
            name:str=None, description:str=None
            ) -> DNASeqRecordWithCoordinates:
        """ concat two sequences into one """
        return self.__class__(
            seq=self.seq + other.seq,
            locations=self.locations + other.locations,
            orf=self.orf,
            id=self.id if id is None else id,
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

    def enzymatic_cleave(self, rule:str, miscleavage:int=2, 
            exception:str=None, start:int=None, end:int=None
            ) -> List[DNASeqRecordWithCoordinates]:
        """ """
        start = 0 if start is None else start
        end = len(self) if end is None else end
        cleave_sites = [start]
        cleave_sites.extend(self.find_all_enzymatic_cleave_sites(rule=rule,
            exception=exception, start=start, end=end))
        cleave_sites.append(end - end % 3 + 1)
        
        fragments = []
        for i in range(1, len(cleave_sites)):
            right_site = cleave_sites[i]
            
            if i == 1:
                left_sites = [0]
            elif i <= miscleavage:
                left_sites = range(i-1, 0, -1)
            else:
                left_sites = range(i-1, i-miscleavage-1, -1)
            
            for left_site in left_sites:
                cdna_seq = self[left_site*3:right_site*3]
                fragments.append(cdna_seq)

        fragments.sort()
        return fragments
    
    def get_query_index(self, ref_index:int) -> int:
        """ Returns the query index wiht a given reference index """
        for location in self.locations:
            if ref_index in location.ref:
                return location.query.start + ref_index - location.ref.start
    
