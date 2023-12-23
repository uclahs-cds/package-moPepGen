""" Module for Amino Acid Record """
from __future__ import annotations
from typing import Iterable, List, Set, Tuple, Union, Pattern
import copy
import re
import regex
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from Bio.Seq import Seq
from moPepGen.SeqFeature import MatchedLocation, FeatureLocation
from moPepGen.aa.expasy_rules import EXPASY_RULES, EXPASY_RULES2, EXPASY_RULES_WINGS_SIZE


class AminoAcidSeqRecord(SeqRecord):
    #pylint: disable=W0223
    """ A AminoAcidSeqRecord holds a protein or peptide sequence
    """
    def __init__(self, seq:Seq, _id:str="<unknown id>",
            name:str="<unknown name>", description:str="<unknown description>",
            gene_id:str=None, transcript_id:str=None, protein_id:str=None,
            gene_name:str=None, **kwargs):
        self.id = None
        self.name = None
        super().__init__(seq, id=_id, name=name, description=description,
            **kwargs)
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.protein_id = protein_id
        self.gene_name = gene_name

    def __getitem__(self, index) -> AminoAcidSeqRecord:
        """"""
        new_seq = self.seq.__getitem__(index)
        new_one = self.__class__(seq=new_seq, _id=self.id, name=self.name,
            description=self.description, gene_id=self.gene_id,
            transcript_id=self.transcript_id, protein_id=self.protein_id)
        return new_one

    def __add__(self, other:AminoAcidSeqRecord) -> AminoAcidSeqRecord:
        """"""
        new_one = super().__add__(other)
        new_one.__class__ = self.__class__
        new_one.gene_id = self.gene_id
        new_one.transcript_id = self.transcript_id
        new_one.protein_id = self.protein_id
        return new_one

    def __hash__(self):
        """ hash """
        return hash(str(self.seq))

    def __eq__(self, other:AminoAcidSeqRecord):
        """ Equal to. It is implemented in this way, because the get_equivalent
        calls this __eq__ instead of the other for reason that I don't
        understand. """
        result = (self.seq == other.seq)
        if result and hasattr(other, 'match'):
            other.match = self
        return result

    def __ne__(self, other:AminoAcidSeqRecord) -> bool:
        """ not equal to """
        return not self == other

    def infer_ids(self, style:str=None) -> str:
        """
        Args:
            source (str): The style of the fasta file header. Either 'GENCODE'
                or 'ENSEMBL'. If None, it will try both. Defaults to None.

        Returns:
            The style ifered.
        """
        if not style:
            try:
                self.infer_ids_gencode()
                return 'GENCODE'
            except ValueError:
                pass

            try:
                self.infer_ids_ensembl()
                return 'ENSEMBL'
            except ValueError as e:
                raise ValueError(
                    'Failed to infer gene ID, transcript ID, and protein ID '
                    r'using both ENSEMBL and GENCODE\'s format.'
                ) from e
        elif style == 'GENCODE':
            self.infer_ids_gencode()
        elif style == 'ENSEMBL':
            self.infer_ids_ensembl()
        else:
            raise ValueError(f'style {style} is not supported')
        return style


    def infer_ids_gencode(self):
        """ Infers the gene, transcript, and protein ID from description base
        on the GENCODE format
        """
        p = r'^(?P<protein_id>[A-Z0-9.:_-]+)\|' +\
            r'(?P<transcript_id>[A-Z0-9.:_-]+)\|(?P<gene_id>[A-Z0-9._]+)\|.+'

        match = re.search(p, self.description)
        if not match:
            raise ValueError(
                f"The record doesn't seem to follow the GENCODE format: {self.description}"
            )

        protein_id = match.group('protein_id')
        transcript_id = match.group('transcript_id')
        gene_id = match.group('gene_id')

        self.id = protein_id
        self.name = protein_id
        self.gene_id = gene_id
        self.protein_id = protein_id
        self.transcript_id = transcript_id

    def infer_ids_ensembl(self):
        """ Infers the gene, transcript, and protein ID from description base
        on the ENSEMBL format
        """
        p = r'^(?P<protein_id>.+) pep (\bchromosome\b|\bsupercontig\b|\bscaffold\b):' +\
            r'(?P<genome_build>.+?):(?P<chromosome>.+?):(?P<position>\d+:\d+):(?P<strand>.+) ' +\
            r'gene:(?P<gene_id>.+) transcript:(?P<transcript_id>\S+) .+'

        match = re.search(p, self.description)
        if not match:
            raise ValueError(f"doesn't seem to follow the ENSEMBL format: {self.description}")

        protein_id = self.id
        gene_id = match.group('gene_id')
        transcript_id = match.group('transcript_id')

        self.gene_id = gene_id.split('.')[0]
        self.protein_id = protein_id.split('.')[0]
        self.transcript_id = transcript_id.split('.')[0]

    def iter_stop_sites(self) -> Iterable[int]:
        """ Create a generator for stop sites """
        for occ in re.finditer(r'\*', str(self.seq)):
            yield occ.start()

    def iter_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            exception_sites:List[int]=None) -> Iterable[int]:
        """ Create a generator of the cleave sites """
        rule = EXPASY_RULES[rule]
        exception = EXPASY_RULES.get(exception, exception)
        if exception_sites is None:
            exception_sites = [] if exception is None else \
                [x.end() for x in re.finditer(exception, str(self.seq))]

        for it in re.finditer(rule, str(self.seq)):
            s = it.end()
            if s not in exception_sites:
                yield it.end()

    def iter_enzymatic_cleave_sites_with_range(self, rule:str, exception:str=None,
            exception_sites:List[int]=None
            ) -> Iterable[Tuple[int, Union[Tuple[int,int], None]]]:
        """ Create a generator of the cleave sites and also yields the full
        range of the pattern recognized during enzymatic cleavage.

        ## Examples:

        In this example, the cleavage site is 5 (after K), and the pattern range
        is (4,6) for KT.
        >>> seq = AminoAcidSeqRecord(Seq('TTTTKTTTT'))
        >>> list(seq.iter_enzymatic_cleave_sites_with_range('trypsin'))
        [(5,(4,6))]

        In this example, the cleavage site is 5 (after R) and the pattern range
        is (3,6) for MRP.
        >>> seq = AminoAcidSeqRecord(Seq('TTTMRPTTT'))
        >>> list(seq.iter_enzymatic_cleave_sites_with_range('trypsin'))
        [(5,(3,6))]
        """
        site_pattern = re.compile(EXPASY_RULES[rule])
        range_pattern = regex.compile(EXPASY_RULES2[rule])
        exception = EXPASY_RULES.get(exception, exception)
        seq = str(self.seq)
        if exception_sites is None:
            exception_sites = [] if exception is None else \
                [x.end() for x in re.finditer(exception, seq)]

        sites = [x.end() for x in site_pattern.finditer(seq)]
        ranges = [(x.start(), x.end()) for x in range_pattern.finditer(seq, overlapped=True)]

        if len(sites) != len(ranges):
            raise ValueError(
                f"Inconsistent cleavage sites found. seq={seq}, "
                f"sites={sites}, ranges={ranges}"
            )
        for s, r in zip(sites, ranges):
            if s not in exception_sites:
                yield s, r

    def iter_enzymatic_cleave_sites_with_range_local(self, rule:str, exception:str=None,
            exception_sites:List[int]=None
            ) -> Iterable[Tuple[int, Union[Tuple[int,int], None]]]:
        """ Iter enzymatic cleavage sites with pattern matching ranges using the
        local matching method. """
        rule_pattern = re.compile(EXPASY_RULES[rule])
        exception = EXPASY_RULES.get(exception, exception)
        wings_size = EXPASY_RULES_WINGS_SIZE[rule]
        seq = str(self.seq)
        seq_len = len(seq)
        if exception_sites is None:
            exception_sites = [] if exception is None else \
                [x.end() for x in re.finditer(exception, seq)]

        sites = [x.end() for x in rule_pattern.finditer(seq)]

        if wings_size[0] <= 0 and wings_size[1] <= 0:
            raise ValueError(
                "Invalid enzyme pattern with the size being 0 for both wings:" +
                f"{rule}: {EXPASY_RULES[rule]}"
            )

        for s in sites:
            if s in exception_sites:
                continue
            r = self.get_local_matched_range(
                seq=seq, site=s, p=rule_pattern, wings_size=wings_size, seq_len=seq_len
            )
            yield s, r

    @staticmethod
    def get_local_matched_range(seq:str, site:int, p:Pattern[str],
            wings_size:Tuple[int, int], seq_len:int=None):
        """ Get the local matched pattern range from a sequence
        with a given cleavage site. """
        if not seq_len:
            seq_len = len(seq)

        upper = max(site - wings_size[0], 0)
        lower = min(site + wings_size[1], len(seq))

        if wings_size[0] >= wings_size[1]:
            ucur = site - 1
            lcur = site
        else:
            ucur = site
            lcur = site + 1

        pattern_found: bool = False
        while True:
            if not upper <= ucur < lcur <= lower:
                break
            local_seg = seq[ucur:lcur]
            if p.search(local_seg):
                pattern_found = True
                break
            if (site - ucur) <= (lcur - site) or lcur == lower:
                ucur -= 1
            else:
                lcur += 1
        if not pattern_found:
            raise ValueError(f"Cannot extract matched pattern at position {site} from {seq}")
        return (ucur, lcur)

    def get_enzymatic_cleave_exception_sites(self, exception:str) -> Iterable[int]:
        """ Find all enzymatic cleavage exception sites """
        return [] if exception is None else \
            [x.end() for x in re.finditer(exception, str(self.seq))]

    def find_first_cleave_or_stop_site(self, rule:str, exception:str=None,
            exception_sites:List[int]=None) -> Iterable[int]:
        """ Create a generator for cleave or stop sites """
        it_cleavage = self.iter_enzymatic_cleave_sites(
            rule=rule, exception=exception, exception_sites=exception_sites
        )
        sites = []
        first_cleave_site = next(it_cleavage, None)
        if first_cleave_site is not None:
            sites.append(first_cleave_site)
        it_stop = self.iter_stop_sites()
        first_stop_site = next(it_stop, None)
        if first_stop_site is not None:
            if first_stop_site == 0:
                if len(self.seq) == 1:
                    return -1
                sites.append(first_stop_site + 1)
            else:
                sites.append(first_stop_site)
        if not sites:
            return -1
        return min(sites)

    def find_first_cleave_or_stop_site_with_range(self, rule:str, exception:str=None,
            exception_sites:List[int]=None
            ) -> Tuple[int, Union[Tuple[int,int], None]]:
        """ Create a generator for cleave or stop sites """
        it_cleavage = self.iter_enzymatic_cleave_sites_with_range(
            rule=rule, exception=exception, exception_sites=exception_sites
        )
        sites:List[Tuple[int, Union[Tuple[int,int], None]]] = []
        first_cleave_site, first_range = next(it_cleavage, (None, None))
        if first_cleave_site is not None:
            sites.append((first_cleave_site, first_range))
        it_stop = self.iter_stop_sites()
        first_stop_site = next(it_stop, None)
        if first_stop_site is not None:
            if first_stop_site == 0:
                if len(self.seq) == 1:
                    return -1, None
                sites.append((first_stop_site + 1, None))
            else:
                sites.append((first_stop_site, None))
        if not sites:
            return -1, None
        return min(sites, key=lambda x: x[0])

    def find_first_enzymatic_cleave_site(self, rule:str, exception:str=None,
            start:int=0) -> int:
        """ Find the first enzymatic cleave site """
        it = self[start:]\
            .iter_enzymatic_cleave_sites(rule=rule, exception=exception)
        try:
            return next(it) + start
        except StopIteration:
            return -1

    def find_all_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            ) -> List[int]:
        """ Find all enzymatic cleave sites. """
        return list(self.iter_enzymatic_cleave_sites(rule=rule,
            exception=exception))

    def find_all_enzymatic_cleave_sites_with_ranges(self, rule:str, exception:str=None,
            ) -> List[Tuple[int, Tuple[int, int]]]:
        """ Find all enzymatic cleave sites. """
        return list(self.iter_enzymatic_cleave_sites_with_range(rule=rule,
            exception=exception))

    def find_all_start_sites(self) -> List[int]:
        """ Find all start positions """
        return [x.start() for x in re.finditer('M', str(self.seq))]

    def find_all_cleave_and_stop_sites(self, rule:str, exception:str=None,
            exception_sites:List[int]=None) -> List[int]:
        """ Find all enzymatic lceave sites and stop sites """
        cleavage_sites = list(self.iter_enzymatic_cleave_sites(rule=rule,
            exception=exception, exception_sites=exception_sites))
        stop_sites_start = list(self.iter_stop_sites())
        stop_sites_end = [i + 1 for i in stop_sites_start if i < len(self.seq) - 1]
        stop_sites_start = [i for i in stop_sites_start if i > 0]
        sites = cleavage_sites + stop_sites_start + stop_sites_end
        sites = list(set(sites))
        sites.sort()
        return sites

    def find_all_cleave_and_stop_sites_with_range(self, rule:str,
            exception:str=None, exception_sites:List[int]=None
            ) -> List[Tuple[int,Union[Tuple[int,int], None]]]:
        """ Find all enzymatic cleavage sites with ranges and stop sites """
        sites = list(
            self.iter_enzymatic_cleave_sites_with_range(
                rule=rule, exception=exception,
                exception_sites=exception_sites
            )
        )
        sites_mapper = dict(sites)
        stop_sites_start = list(self.iter_stop_sites())
        stop_sites_end = [i + 1 for i in stop_sites_start if i < len(self.seq) - 1]
        stop_sites_start = [i for i in stop_sites_start if i > 0]

        for s in stop_sites_start:
            sites_mapper[s] = (s, s+1)
        for s in stop_sites_end:
            sites_mapper[s] = (s-1, s)

        sites = list(sites_mapper.items())
        sites.sort(key=lambda x: x[0])
        return sites

    def enzymatic_cleave(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500.0, min_length:int=7,
            max_length:int=25, cds_start_nf:bool=False
            )->Set[AminoAcidSeqRecord]:
        """ Performs enzymatic cleave """
        peptides = []
        sites = [0]
        sites += self.find_all_enzymatic_cleave_sites(
            rule=rule, exception=exception)
        sites.append(len(self))
        start = 0

        def update_peptides(peptide):
            if 'X' in peptide.seq:
                return
            mol_wt = SeqUtils.molecular_weight(peptide.seq, 'protein')
            weight_flag = mol_wt > min_mw
            length_flag = len(peptide.seq) >= min_length \
                and len(peptide.seq) <= max_length
            if weight_flag and length_flag:
                peptides.append(peptide)

        while start < len(sites) - 1:
            end = start + 1
            while end - start - 1 <= miscleavage and end < len(sites):
                peptide = self[sites[start]:sites[end]]
                if start == 0 and not cds_start_nf and peptide.seq.startswith('M'):
                    update_peptides(peptide[1:])
                update_peptides(peptide)
                end += 1
            start += 1
        return peptides


class AminoAcidSeqRecordWithCoordinates(AminoAcidSeqRecord):
    """ Amino acid sequence record with coordinates """
    def __init__(self, seq:Seq, *args,
            locations:List[MatchedLocation]=None, orf:FeatureLocation=None, **kwargs ):
        """ Constract a DNASeqRecordWithCoordinates object.

        Args:
            seq (Seq): The amino acid sequence.
            locations (List[MatchedLocation]): The locations where this
                sequence is aligned tp.
            orf (FeatureLocation): The open reading frame start and end.
        """
        super().__init__(seq=seq, *args, **kwargs)
        self.locations = locations or []
        # query index
        self.orf = orf

    def __add__(self, other:AminoAcidSeqRecordWithCoordinates
            ) -> AminoAcidSeqRecordWithCoordinates:
        """ add operation """
        new = super().__add__(other)
        new.__class__ = AminoAcidSeqRecordWithCoordinates
        left_locs = copy.copy(self.locations)
        right_locs = copy.copy(other.locations)
        if left_locs and right_locs:
            lhs = self.locations[-1]
            rhs = other.locations[0].shift(len(self))
            if lhs.query.end == rhs.query.start \
                    and lhs.ref.end == rhs.ref.start \
                    and lhs.get_ref_dna_end() == rhs.get_ref_dna_start():
                query = FeatureLocation(
                    start=lhs.query.start, end=rhs.query.end,
                    reading_frame_index=lhs.query.reading_frame_index,
                    start_offset=lhs.query.start_offset,
                    end_offset=rhs.query.end_offset
                )
                ref = FeatureLocation(
                    start=lhs.ref.start, end=rhs.ref.end, seqname=lhs.ref.seqname,
                    reading_frame_index=lhs.ref.reading_frame_index,
                    start_offset=lhs.ref.start_offset,
                    end_offset=rhs.ref.end_offset
                )
                new_loc = MatchedLocation(query=query, ref=ref)
                right_locs.pop(0)
                left_locs[-1] = new_loc

        new.locations = left_locs + [loc.shift(len(self)) for loc in right_locs]
        new.orf = self.orf
        return new

    def __getitem__(self, index)->AminoAcidSeqRecordWithCoordinates:
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
            _id=self.id,
            name=self.name,
            description=self.description
        )

    def __hash__(self):
        """hash"""
        return hash(self.seq)

    def get_query_index(self, ref_index:int, seqname:str, reading_frame:int=None) -> int:
        """ Returns the query index wiht a given reference index """
        for location in self.locations:
            if location.ref.seqname != seqname:
                continue
            if reading_frame is not None and location.query.reading_frame_index != reading_frame:
                continue
            if ref_index in location.ref:
                return location.query.start + ref_index - location.ref.start
        return -1
