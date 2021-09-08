""" module for peptide pool splitter """
from __future__ import annotations
from typing import Dict, IO, Iterable, List, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen.seqvar import GVFMetadata
from moPepGen import err, seqvar, VARIANT_PEPTIDE_SOURCE_DELIMITER, \
    SPLIT_DATABASE_KEY_SEPARATER
from .VariantPeptidePool import VariantPeptidePool


if TYPE_CHECKING:
    from .AminoAcidSeqRecord import AminoAcidSeqRecord

NONCODING_SOURCE = 'Noncoding'

class VariantSourceSet(set):
    """ Variant source set. This is a class of ordered set.

    Example:
        >>> VariantSourceSet.set_levels({'A':0, 'B':1, 'C':2})
        >>> print(VariantSourceSet(['A']) < VariantSourceSet(['B']))
        True

        >>> print(VariantSourceSet(['A', 'B']) > VariantSourceSet(['B']))
        True
    """
    levels_map = {}
    levels = []

    def __init__(self, *args):
        """ constructor """
        for it in list(*args):
            self.validate(it)
        super().__init__(*args)

    def validate(self, it):
        """ Validate element """
        if it not in self.levels_map:
            raise ValueError(f'No defined level found for {it}')

    def add(self, element):
        """ override the add method """
        self.validate(element)
        super().add(element)

    @classmethod
    def set_levels(cls, levels:Dict[str,int]):
        """ set levels. This is to make sure that all instances share the same
        order. """
        cls.levels_map = levels
        cls.levels = [y[0] for y in sorted(cls.levels_map.items(), key=lambda x:x[1])]

    @classmethod
    def reset_levels(cls):
        """ reset levels """
        cls.levels_map = {}
        cls.levels = []

    def __str__(self) -> str:
        """ str """
        sorted_list = [x for x in self.levels if x in self]
        return SPLIT_DATABASE_KEY_SEPARATER.join(sorted_list)

    def __gt__(self, other:VariantSourceSet) -> bool:
        """ greater than """
        if self == other:
            return False
        if len(self) > len(other):
            return True
        if len(self) < len(other):
            return False
        this = self.to_int()
        that = other.to_int()
        for i, j in zip(this, that):
            if i > j:
                return True
            if i < j:
                return False
        return False

    def __ge__(self, other:VariantSourceSet) -> bool:
        """ greater or equal to """
        return self == other or self > other

    def __lt__(self, other:VariantSourceSet) -> bool:
        """ less then """
        return not self >= other

    def __le__(self, other:VariantSourceSet) -> bool:
        """ less or equal to """
        return not self > other

    def to_int(self, sort=True) -> Iterable[int]:
        """ to int """
        source_int = {self.levels_map[x] for x in self}
        if sort:
            source_int = list(source_int)
            source_int.sort()
        return source_int

class VariantPeptideInfo():
    """ Variant peptide label. This is a helper class in order to sort peptide
    labels easily. """
    def __init__(self, transcript_id:str, variant_labels:List[str],
            variant_index:int, sources:VariantSourceSet=None):
        """ Constructor """
        self.transcript_id = transcript_id
        self.variant_labels = variant_labels
        self.variant_index = variant_index
        self.sources = sources or VariantSourceSet()

    @staticmethod
    def from_variant_peptide(peptide:AminoAcidSeqRecord, label_map:LabelMap
            ) -> List[VariantPeptideInfo]:
        """ Parse from a variant peptide record """
        info_list = []
        delimiter = VARIANT_PEPTIDE_SOURCE_DELIMITER
        for label in peptide.description.split(delimiter):
            tx_id, *var_ids, var_index = label.split('|')
            info = VariantPeptideInfo(tx_id, var_ids, var_index)
            if not var_ids:
                info.sources.add(NONCODING_SOURCE)
            for var_id in var_ids:
                try:
                    source = label_map[tx_id][var_id]
                except KeyError as e:
                    raise err.VariantSourceNotFoundError(tx_id, var_id) from e
                info.sources.add(source)
            info_list.append(info)
        return info_list

    def __str__(self) -> str:
        """ str """
        x = [self.transcript_id] + self.variant_labels + [self.variant_index]
        return '|'.join(x)

    def __eq__(self, other:VariantPeptideInfo):
        """ equal to """
        return self.transcript_id == other.transcript_id\
            and set(self.variant_labels) == set(other.variant_labels) \
            and self.variant_index == other.variant_index \
            and self.sources == other.sources

    def __ne__(self, other:VariantPeptideInfo):
        """ not equal to """
        return not self == other

    def __gt__(self, other:VariantPeptideInfo):
        """" greater than """
        return self.sources > other.sources

    def __ge__(self, other:VariantPeptideInfo):
        """ greater than or equal to """
        return self > other or self == other

    def __lt__(self, other:VariantPeptideInfo):
        """ less than """
        return not self >=  other

    def __le__(self, other:VariantPeptideInfo):
        """ less than or equal to """
        return not self > other

LabelMap = Dict[str,Dict[str,str]]
Databases = Dict[str,VariantPeptidePool]

class PeptidePoolSplitter():
    """ PeptidePoolSplitter to split peptide pool into separate databases.

    Attributes:
        peptides (VariantPeptidePool): The variant peptide pool.
        order (Dict[str,str]): Order of sources.
        label_map (Dict[str,Dict[str,str]]): This a nested dict of transcripts,
            variant lables, and their corresponding source type. This is to
            quickly retrieve the source type of a given variant label.
        group_map (Dict[str,str]): A map of certain variant sources being
            grouped together.
        databases (Dict[str,VariantPeptidePool]): This is the splitted
            databases. It is a dict that keys are the string representation of
            variant sources, and values are the splitted VariantPeptidePool.
        sources (Set[str]): Variant sources.
    """
    def __init__(self, peptides:VariantPeptidePool=None, order:Dict[str,int]=None,
            label_map:LabelMap=None, group_map:Dict[str,str]=None,
            databases:Databases=None, sources:Set[str]=None):
        self.peptides = peptides
        self.databases = databases or {}
        self.label_map = label_map or {}
        self.group_map = group_map or {}
        self.order = order or {}
        self.sources = sources or set()

    def append_order_noncoding(self):
        """ Add noncoding to the end of order """
        source = NONCODING_SOURCE
        self.sources.add(source)
        if source not in self.order:
            self.append_order(source)

    def append_order(self, source:str):
        """ Add a group to the end of the order. """
        if source in self.order:
            raise ValueError(f'The source {source} already has an order.')
        self.order[source] = max(self.order.values()) + 1 if self.order else 0

    def load_database(self, handle:IO) -> None:
        """ load peptide database """
        pool = VariantPeptidePool.load(handle)
        if not self.peptides:
            self.peptides = pool
            return
        for peptide in pool.peptides:
            self.peptides.add_peptide(peptide, None, skip_checking=True)

    def load_database_noncoding(self, handle:IO) -> None:
        """ Load noncoding peptide database """
        self.load_database(handle)
        self.append_order_noncoding()

    def load_gvf(self, handle:IO) -> None:
        """ Load variant lables from GVF file. """
        metadata = GVFMetadata.parse(handle)
        variants = seqvar.io.parse(handle)
        source = metadata.source

        if source in self.group_map:
            source = self.group_map[source]

        self.sources.add(source)

        if source not in self.order:
            self.append_order(source)

        for variant in variants:
            tx_id = variant.attrs['TRANSCRIPT_ID']
            if tx_id not in self.label_map:
                self.label_map[tx_id] = {}
            if variant.id not in self.label_map[tx_id]:
                self.label_map[tx_id][variant.id] = source

    def add_peptide_to_database(self, database_key:str,
            peptide:AminoAcidSeqRecord) -> None:
        """ Add peptide to database """
        if database_key not in self.databases:
            self.databases[database_key] = VariantPeptidePool()
        self.databases[database_key].peptides.add(peptide)

    def split(self, max_groups:int, additional_split:List[Set]):
        """ Split peptide pool into separate databases """
        VariantSourceSet.set_levels(self.order)
        delimiter = VARIANT_PEPTIDE_SOURCE_DELIMITER
        additional_split = [VariantSourceSet(x) for x in additional_split]
        for peptide in self.peptides.peptides:
            peptide_infos = VariantPeptideInfo.from_variant_peptide(
                peptide, self.label_map)
            peptide_infos.sort()

            peptide.description = delimiter.join([str(x) for x in peptide_infos])
            peptide.id = peptide.description
            peptide.name = peptide.description

            if len(peptide_infos[0].sources) <= max_groups:
                database_key = str(peptide_infos[0].sources)
                self.add_peptide_to_database(database_key, peptide)
            else:
                has_additional_splitting = False
                for additional_set in additional_split:
                    ind = get_index_contains_source_set(peptide_infos, additional_set)
                    if ind is not None:
                        info = peptide_infos.pop(ind)
                        peptide_infos.insert(0, info)
                        label = delimiter.join([str(x) for x in peptide_infos])
                        peptide.description = label
                        peptide.id = peptide.description
                        peptide.name = peptide.description

                        database_key = self.get_additional_database_key(additional_set)
                        self.add_peptide_to_database(database_key, peptide)
                        has_additional_splitting = True
                        break
                if not has_additional_splitting:
                    database_key = self.get_remaining_database_key()
                    self.add_peptide_to_database(database_key, peptide)

    @staticmethod
    def get_additional_database_key(additional:str) -> str:
        """ Get the database key for priority """
        return f'{str(additional)}{SPLIT_DATABASE_KEY_SEPARATER}additional'

    @staticmethod
    def get_remaining_database_key() -> str:
        """ Get the database key for remanent peptides """
        return 'Remaining'

    @staticmethod
    def get_database_filename(database_key:str, prefix:str) -> str:
        """ Get the database filename """
        return f"{prefix}_{database_key}.fasta"

    def write(self, output_prefix:str) -> None:
        """ Write databases to disk """
        output_dir = Path(output_prefix).parent
        output_dir.mkdir(exist_ok=True)
        filename_prefix = str(Path(output_prefix).name)
        for database_key, database in self.databases.items():
            filename = self.get_database_filename(database_key, filename_prefix)
            database.write(output_dir/filename)


def get_index_contains_source_set(peptide_infos:List[VariantPeptideInfo],
        target_set:VariantSourceSet) -> List[int]:
    """ Find the indices of a list of source sets,  """
    for i, info in enumerate(peptide_infos):
        if target_set.issubset(info.sources):
            return i
    return None
