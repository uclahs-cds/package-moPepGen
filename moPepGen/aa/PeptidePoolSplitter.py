""" module for peptide pool splitter """
from __future__ import annotations
from typing import Dict, Iterable, List, Set, TYPE_CHECKING, Union
from moPepGen.seqvar import GVFMetadata
from moPepGen import seqvar, VARIANT_PEPTIDE_DELIMITER
from .VariantPeptidePool import VariantPeptidePool


if TYPE_CHECKING:
    from pathlib import Path
    from .AminoAcidSeqRecord import AminoAcidSeqRecord

NONCODING_SOURCE = 'Noncoding'

class VariantSourceSet(set):
    """ Variant source set """
    levels = []

    @classmethod
    def set_levels(cls, levels:Dict[str,int]):
        """ set levels """
        cls.levels = levels

    def __str__(self) -> str:
        """ str """
        return '+'.join(self)

    def __gt__(self, other:VariantSourceSet) -> bool:
        """ greater than """
        if self == other:
            return False
        if len(self) > len(other):
            return True
        this = self.to_int()
        that = other.to_int()
        for x, y in enumerate(this, that):
            if x > y:
                return True
            if x < y:
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
        for label in peptide.description.split(VARIANT_PEPTIDE_DELIMITER):
            tx_id, *var_ids, _ = label.split('|')
            sources = cls()
            if not len(var_ids):
                sources.add(NONCODING_SOURCE)
            for var_id in var_ids:
                sources.add(label_map[tx_id][var_id])
            source_list.append(sources)
        return source_list

LabelMap = Dict[str,Dict[str,str]]
DatabaseKey = Union[VariantSourceSet, str]
Databases = Dict[DatabaseKey,VariantPeptidePool]

class PeptidePoolSplitter():
    """ PeptidePoolSplitter to split peptide pool into separate databases """
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
        self.order[source] = max(self.order.values()) + 1

    def load_database(self, path:Path) -> None:
        """ load peptide database """
        pool = VariantPeptidePool.load(path)
        if not self.peptides:
            self.peptides = pool
            return
        for peptide in pool.peptides:
            self.peptides.add_peptide(peptide)

    def load_database_noncoding(self, path:Path) -> None:
        """ Load noncoding peptide database """
        self.load_database(path)
        self.append_order_noncoding()

    def load_gvf(self, path:Path) -> None:
        """ Load variant lables from GVF file. """
        with open(path, 'rt') as handle:
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
        for group in priority_list:
            if group not in self.sources:
                raise ValueError(f"group {group} not found")
        VariantSourceSet.set_levels(self.order)
        for peptide in self.peptides:
            source_sets = VariantSourceSet.from_variant_peptide(
                peptide, self.label_map)

            order = [i for (v, i) in sorted((v, i) for (i, v) in enumerate(source_sets))]
            source_sets = [source_sets[i] for i in order]

            labels = peptide.description.split(VARIANT_PEPTIDE_DELIMITER)
            labels = [labels[i] for i in order]
            peptide.description = VARIANT_PEPTIDE_DELIMITER.join(labels)
            peptide.id = peptide.description
            peptide.name = peptide.description

            if len(source_sets[0]) < max_groups:
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
                        peptide.description = VARIANT_PEPTIDE_DELIMITER.join(labels)
                        peptide.id = peptide.description
                        peptide.name = peptide.description

                        database_key = self.get_priority_database_key(priority)
                        self.add_peptide_to_database(database_key, peptide)
                        has_priority = True
                        break
                if not has_priority:
                    database_key = self.get_remanent_database_key()
                    self.add_peptide_to_database(database_key, peptide)

    @staticmethod
    def get_priority_database_key(priority:str) -> str:
        """ Get the database key for priority """
        return f'{priority}+'

    @staticmethod
    def get_remanent_database_key() -> str:
        """ Get the database key for remanent peptides """
        return 'Remanent'

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
