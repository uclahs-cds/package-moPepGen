""" module for peptide pool summarizer """
from __future__ import annotations
import itertools
import statistics
from typing import Dict, IO, List, Set, FrozenSet, Tuple, Optional
import matplotlib.pyplot as plt
from moPepGen import seqvar, constant
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from moPepGen.aa.VariantPeptideLabel import VariantPeptideInfo, \
    VariantSourceSet, LabelSourceMapping
from moPepGen.aa.VariantPeptidePool import VariantPeptidePool
from moPepGen.seqvar.GVFMetadata import GVFMetadata


SOURCES_INTERNAL = [
    constant.SOURCE_NOVEL_ORF,
    constant.SOURCE_SEC_TERMINATION,
    constant.SOURCE_CODON_REASSIGNMENT
]

MUTUALLY_EXCLUSIVE_PARSERS:Dict[str,List[str]] = {
    'parseSTARFusion': [
        'parseFusionCatcher',
        'parseArriba',
        'parseCIRCexplorer',
        'parseRMATS'
    ],
    'parseFusionCatcher': [
        'parseSTARFusion',
        'parseArriba',
        'parseCIRCexplorer',
        'parseRMATS'
    ],
    'parseArriba': [
        'parseSTARFusion',
        'parseFusionCatcher',
        'parseCIRCexplorer',
        'parseRMATS'
    ],
    'parseCIRCexplorer': [
        'parseSTARFusion',
        'parseFusionCatcher',
        'parseArriba',
        'parseRMATS'
    ],
    'parseRMATS': [
        'parseSTARFusion',
        'parseFusionCatcher',
        'parseArriba',
        'parseCIRCexplorer'
    ]
}

# BPF default color scheme.
# Hopefully no one is using more than 11 miscleavages
COLOR_SCHEME = [
    "#FFA500", "#458B00", "#68228B", "#FFD700", "#1E90FF", "#CD2626", "#9ACD32"
    "#FF7F00", "#473C8B", "#43CD80", "#CD3278", "#00C5CD"
]

class NoncanonicalPeptideSummaryTable():
    """ Summary table for noncaonincal peptides

    Attributes:
        data (Dict[FrozonSet[str], Dict[str,int]]): Mapping from variant source
            to number of miscleavages to count of peptides.
        max_misc (int): Max number of miscleavages allowed.
    """
    KEY_N_TOTAL = 'n_total'
    def __init__(self, data:Dict[FrozenSet[str], Dict[str,int]]=None,
            sources:Set[str]=None, max_misc:int=0):
        """ constructor """
        # mapping: source -> misc -> count
        self.data = data or {}
        self.sources = sources or set()
        self.max_misc = max_misc

    @staticmethod
    def get_key_x_misc(x:int):
        """ get the table column name for number of peptides per combination of
        sources with x number of miscleavages """
        return f"n_{x}_misc"

    def get_n_total(self, key:str) -> int:
        """ gety the table column name for total number of peptides per
        combination of sources """
        return self.data[key][self.KEY_N_TOTAL]

    def get_n_x_misc(self, key:str, misc:int) -> int:
        """ get the number of peptides with x number of miscalvages from the
        summary table """
        entry = self.data[key]
        key = self.get_key_x_misc(misc)
        return entry.get(key, 0)

    def add_entry(self,
            seq:AminoAcidSeqRecord,
            label_map:LabelSourceMapping,
            group_map:Dict[str,str],
            tx2gene:Dict[str,str],
            coding_tx:Set[str],
            enzyme:str) -> None:
        """ Add peptide entry to the summary table """
        peptide_labels = VariantPeptideInfo.from_variant_peptide(
            peptide=seq, tx2gene=tx2gene, coding_tx=coding_tx,
            label_map=label_map, group_map=group_map,
            check_source=True
        )
        peptide_labels.sort()
        sources = frozenset(peptide_labels[0].sources)
        self.increment_total(sources)

        exception = 'trypsin_exception' if enzyme == 'trypsin' else None
        misc = len(seq.find_all_enzymatic_cleave_sites(enzyme, exception))
        self.increment_misc(sources, misc)

    def add_source(self, source:str) -> None:
        """ add source """
        self.sources.add(source)

    def _increment(self, sources:FrozenSet[str], key):
        """ Increment """
        for source in sources:
            self.add_source(source)

        if sources not in self.data:
            self.data[sources] = {key: 1}
        elif key not in self.data[sources]:
            self.data[sources][key] = 1
        else:
            self.data[sources][key] += 1

    def increment_total(self, sources:FrozenSet[str]):
        """ Increment the total number of peptides for the given combination of
        sources """
        self._increment(sources, self.KEY_N_TOTAL)

    def increment_misc(self, sources:FrozenSet[str], misc:int):
        """ Increment the number of peptides with x miscleavages for the given
        combination of sources """
        self.max_misc = max(self.max_misc, misc)
        self._increment(sources, self.get_key_x_misc(misc))

    def get_keys(self) -> Tuple[str]:
        """ Get all column names """
        return ('n_total', *[self.get_key_x_misc(i) for i in range(self.max_misc + 1)])

    def get_stringified_summary_entry(self, key:str, order:Dict[str,int],
            sep:str='\t') -> str:
        """ Get a row from the summary table for a given combination of sources. """
        rowname = '-'.join(sorted(key, key=lambda x: order[x]))
        if key not in self.data:
            entry = [rowname, '0']
            for _ in range(self.max_misc + 1):
                entry.append('0')
        else:
            entry = [rowname, str(self.get_n_total(key))]
            for x in range(self.max_misc + 1):
                entry.append(str(self.get_n_x_misc(key, x)))
        return sep.join(entry)

    def has_source(self, source:str) -> bool:
        """ Checks if the given source exist in the summary table """
        return source in self.sources

class PeptidePoolSummarizer():
    """ Summarize the variant peptide pool called by moPepGen callVariant. """
    def __init__(self,
            peptides:VariantPeptidePool=None,
            summary_table:NoncanonicalPeptideSummaryTable=None,
            label_map:LabelSourceMapping=None,
            order:Dict[str,int]=None,
            group_map:Dict[str,str]=None,
            source_parser_map:Dict[str,str]=None,
            ignore_missing_source:bool=False):
        """ """
        self.peptides = peptides
        self.summary_table = summary_table or NoncanonicalPeptideSummaryTable()
        self.label_map = label_map or LabelSourceMapping()
        self.order = order or {}
        self.group_map = group_map or {}
        self.source_parser_map = source_parser_map or {}
        self.ignore_missing_source = ignore_missing_source

    def get_reversed_group_map(self) -> Dict[str, List[str]]:
        """ Reverse group map """
        group_map:Dict[str, List[str]] = {}
        for k,v in self.group_map.items():
            if v in group_map:
                group_map[v].append(k)
            else:
                group_map[v] = [k]
        return group_map

    def get_parsers_from_source(self, source:str) -> Set[str]:
        """ Get parsers from source """
        group_map = self.get_reversed_group_map()
        parsers = set()
        # pylint: disable=R1715
        if source in group_map:
            sources = group_map[source]
        else:
            sources = [source]

        for s in sources:
            if s in SOURCES_INTERNAL:
                continue
            parser = self.source_parser_map[s]
            parsers.add(parser)
        return parsers

    def append_order_internal_sources(self):
        """ Add internal sources that are not present in any GTFs, including
        novel ORF, sec termination, and codon reassignment. """
        for source in SOURCES_INTERNAL:
            if source in self.group_map:
                source = self.group_map[source]
            if source not in self.order:
                self.append_order(source)

    def append_order(self, source:str):
        """ Add a group to the end of the order. """
        if source in self.order:
            return
        if source in self.group_map:
            source = self.group_map[source]
        if source in self.order:
            return
        self.order[source] = max(self.order.values()) + 1 if self.order else 0

    def load_database(self, handle:IO) -> None:
        """ load peptide database """
        pool = VariantPeptidePool.load(handle)
        if not self.peptides:
            self.peptides = pool
            return
        for peptide in pool.peptides:
            self.peptides.add_peptide(peptide, None, skip_checking=True)

    def update_label_map(self, handle:IO):
        """ Read GVF files and add to the label source map. """
        metadata = GVFMetadata.parse(handle)
        self.append_order(metadata.source)
        self.source_parser_map[metadata.source] = metadata.parser
        for gene_id, _, label in seqvar.io.parse_label(handle):
            self.label_map.add_record(gene_id, label, metadata.source)

    def count_peptide_source(self, tx2gene:Dict[str,str], coding_tx:Set[str], enzyme:str):
        """ Count number of peptides in each source or combinations of sources. """
        VariantSourceSet.set_levels(self.order)
        for seq in self.peptides.peptides:
            self.summary_table.add_entry(
                seq=seq,
                label_map=self.label_map,
                group_map=self.group_map,
                tx2gene=tx2gene,
                coding_tx=coding_tx,
                enzyme=enzyme
            )

    def contains_exclusive_sources(self, sources:Set[str]) -> bool:
        """ Checks whether the given sources contains any mutually exclusive
        parsers. """
        group_map = self.get_reversed_group_map()
        ungrouped_sources = set()
        for source in sources:
            if source in group_map:
                ungrouped_sources.update(group_map[source])
            else:
                ungrouped_sources.add(source)

        # Two source groups are exclusive only when every source in a group A
        # has at least one mutually exclusive group in group B, and also the
        # other way around.

        for source in sources:
            parsers = self.get_parsers_from_source(source)
            if not parsers:
                continue

            for other in sources:
                other_parsers = self.get_parsers_from_source(other)
                if not other_parsers:
                    continue
                is_exclusive = True
                for x in parsers:
                    if x not in MUTUALLY_EXCLUSIVE_PARSERS:
                        is_exclusive = False
                        break
                    for y in other_parsers:
                        if y not in MUTUALLY_EXCLUSIVE_PARSERS:
                            is_exclusive = False
                            break
                        if x not in MUTUALLY_EXCLUSIVE_PARSERS[y]:
                            is_exclusive = False
                            break

                if is_exclusive:
                    return True

        return False

    def write_summary_table(self, handle:IO):
        """ Write summary table to output. """
        summary_keys = self.summary_table.get_keys()
        header = '\t'.join(['sources', *summary_keys])
        handle.write(header + '\n')
        sources = [it[0] for it in sorted(self.order.items(), key=lambda x:x[1])
                    if not isinstance(it[0], frozenset)]
        for i in range(len(sources)):
            for comb in itertools.combinations(sources, i + 1):
                if self.ignore_missing_source:
                    # Ignore it if the source isn't present in any GVF given.
                    if any(not self.summary_table.has_source(k) for k in comb):
                        continue
                if self.contains_exclusive_sources(comb):
                    continue
                key = frozenset(comb)
                record = self.summary_table.get_stringified_summary_entry(key, self.order)
                handle.write(record + '\n')

    def create_barplot(self, width:float=6, height:float=8, ax:Optional[plt.Axes]=None,
            scale:str=None) -> plt.Axes:
        """ Make a barplot of the summarized data. """
        # misc -> source -> count
        data:Dict[int,List[int]] = {}
        keys:List[str] = []
        sources = [it[0] for it in sorted(self.order.items(), key=lambda x:x[1])]
        for i in range(len(sources)):
            for comb in itertools.combinations(sources, i + 1):
                if self.ignore_missing_source:
                    # Ignore it if the source isn't present in any GVF given.
                    if any(not self.summary_table.has_source(k) for k in comb):
                        continue
                if self.contains_exclusive_sources(comb):
                    continue
                key = frozenset(comb)
                if key not in self.summary_table.data:
                    continue

                keys.append('-'.join(key))
                for x in range(self.summary_table.max_misc + 1):
                    n = self.summary_table.get_n_x_misc(key, x)
                    if x not in data:
                        data[x] = []
                    data[x].append(n)

        totals = []
        for x in data.values():
            if not totals:
                totals = x
            else:
                totals = [i+j for i,j in zip(totals, x)]
                mean = statistics.mean(totals)
                median = statistics.median(totals)

        # Unless specified, log scale is used if the summary distribution is skewed.
        if scale not in ['normal', 'log']:
            scale = 'normal' if mean / median <= 2.5 else 'log'

        if ax is None:
            _, ax = plt.subplots(figsize=(width, height))

        if scale == 'log':
            ax.barh(
                y=list(reversed(keys)), width=list(reversed(totals))
            )
            ax.set_xscale('log')
        else:
            offset = None
            for i,n in enumerate(data.keys()):
                vals = list(reversed(data[n]))
                if offset:
                    vals = [x + y for x,y in zip(vals, offset)]
                ax.barh(
                    y=list(reversed(keys)), width=vals, left=offset, label=n,
                    color=COLOR_SCHEME[i]
                )
                offset = vals
                ax.legend(title='Miscleavages')
        ax.set_xlabel('Number of Peptides')
        return ax
