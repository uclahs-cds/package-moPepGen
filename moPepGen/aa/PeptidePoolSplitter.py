""" module for peptide pool splitter """
from __future__ import annotations
from typing import Dict, IO, List, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen.seqvar import GVFMetadata
from moPepGen import seqvar, circ, VARIANT_PEPTIDE_SOURCE_DELIMITER, \
    SPLIT_DATABASE_KEY_SEPARATER
from .VariantPeptidePool import VariantPeptidePool
from .VariantPeptideLabel import VariantPeptideInfo, VariantSourceSet, \
    LabelSourceMapping

if TYPE_CHECKING:
    from .AminoAcidSeqRecord import AminoAcidSeqRecord
    from moPepGen.gtf import GenomicAnnotation

NONCODING_SOURCE = 'Noncoding'

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
            label_map:LabelSourceMapping=None, group_map:Dict[str,str]=None,
            databases:Databases=None, sources:Set[str]=None):
        self.peptides = peptides
        self.databases = databases or {}
        self.label_map = label_map or LabelSourceMapping()
        self.group_map = group_map or {}
        self.order = order or {}
        self.sources = sources or set()

    def append_order_noncoding(self):
        """ Add noncoding to the end of order """
        source = NONCODING_SOURCE
        if source in self.sources:
            return
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

    def load_gvf(self, handle:IO) -> None:
        """ Load variant lables from GVF file. """
        metadata = GVFMetadata.parse(handle)
        source = metadata.source

        if source in self.group_map:
            source = self.group_map[source]

        self.sources.add(source)

        if source not in self.order:
            self.append_order(source)

        if metadata.parser == 'parseCIRCexplorer':
            for record in circ.io.parse(handle):
                self.label_map.add_circ_rna(record, source)
        else:
            for variant in seqvar.io.parse(handle):
                self.label_map.add_variant(variant, source)

    def add_peptide_to_database(self, database_key:str,
            peptide:AminoAcidSeqRecord) -> None:
        """ Add peptide to database """
        if database_key not in self.databases:
            self.databases[database_key] = VariantPeptidePool()
        self.databases[database_key].peptides.add(peptide)

    def split(self, max_groups:int, additional_split:List[Set],
            anno:GenomicAnnotation):
        """ Split peptide pool into separate databases """
        self.append_order_noncoding()
        VariantSourceSet.set_levels(self.order)
        delimiter = VARIANT_PEPTIDE_SOURCE_DELIMITER
        additional_split = [VariantSourceSet(x) for x in additional_split]
        for peptide in self.peptides.peptides:
            peptide_infos = VariantPeptideInfo.from_variant_peptide(
                peptide, anno, self.label_map)
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
