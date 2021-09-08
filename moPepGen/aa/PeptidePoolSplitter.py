""" module for peptide pool splitter """
from __future__ import annotations
from typing import Dict, IO, Iterable, List, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen.seqvar import GVFMetadata
from moPepGen import seqvar, VARIANT_PEPTIDE_SOURCE_DELIMITER
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
    levels = {}

    def __init__(self, *args):
        """ constructor """
        for it in list(*args):
            self.validate(it)
        super().__init__(*args)

    def validate(self, it):
        """ Validate element """
        if it not in self.levels:
            raise ValueError(f'No defined level found for {it}')

    def add(self, element):
        """ override the add method """
        self.validate(element)
        super().add(element)

    @classmethod
    def set_levels(cls, levels:Dict[str,int]):
        """ set levels. This is to make sure that all instances share the same
        order. """
        cls.levels = levels

    @classmethod
    def reset_levels(cls):
        """ reset levels """
        cls.levels = {}

    def __str__(self) -> str:
        """ str """
        return '+'.join(self)

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
        source_int = {self.levels[x] for x in self}
        if sort:
            source_int = list(source_int)
            source_int.sort()
        return source_int

    @classmethod
    def from_variant_peptide(cls, peptide:AminoAcidSeqRecord,
            label_map:LabelMap) -> List[VariantSourceSet]:
        """ Create VariantSourceSet from AminoAcidSeqRecord """
        source_list = []
        for label in peptide.description.split(VARIANT_PEPTIDE_SOURCE_DELIMITER):
            tx_id, *var_ids, _ = label.split('|')
            sources = cls()
            if not var_ids:
                sources.add(NONCODING_SOURCE)
            for var_id in var_ids:
                sources.add(label_map[tx_id][var_id])
            source_list.append(sources)
        return source_list

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

    def split(self, max_groups:int, priority_list:List[str]):
        """ Split peptide pool into separate databases """
        delimiter = VARIANT_PEPTIDE_SOURCE_DELIMITER
        for group in priority_list:
            if group not in self.sources:
                raise ValueError(f"group {group} not found")
        VariantSourceSet.set_levels(self.order)
        for peptide in self.peptides.peptides:
            source_sets = VariantSourceSet.from_variant_peptide(
                peptide, self.label_map)

            order = [i for _,i in sorted((v,i) for i,v in enumerate(source_sets))]
            source_sets = [source_sets[i] for i in order]

            labels = peptide.description.split(delimiter)
            labels = [labels[i] for i in order]
            peptide.description = delimiter.join(labels)
            peptide.id = peptide.description
            peptide.name = peptide.description

            if len(source_sets[0]) <= max_groups:
                database_key = str(source_sets[0])
                self.add_peptide_to_database(database_key, peptide)
            else:
                has_priority = False
                for priority in priority_list:
                    ind = get_index_with_priority(source_sets, priority)
                    if ind is not None:
                        sources = source_sets.pop(ind)
                        source_sets.insert(0, sources)
                        label = labels.pop(ind)
                        labels.insert(label, 0)
                        peptide.description = delimiter.join(labels)
                        peptide.id = peptide.description
                        peptide.name = peptide.description

                        database_key = self.get_priority_database_key(priority)
                        self.add_peptide_to_database(database_key, peptide)
                        has_priority = True
                        break
                if not has_priority:
                    database_key = self.get_remaining_database_key()
                    self.add_peptide_to_database(database_key, peptide)

    @staticmethod
    def get_priority_database_key(priority:str) -> str:
        """ Get the database key for priority """
        return f'{priority}+'

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


def get_index_with_priority(source_sets:List[VariantSourceSet],
        priority:str) -> List[int]:
    """ Find the indices of sources that contains label from a priority list.
    The list must be already sorted. """
    for i, sources in enumerate(source_sets):
        if priority in sources:
            return i
    return None
