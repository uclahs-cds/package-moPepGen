""" module for peptide pool summarizer """
from __future__ import annotations
import itertools
from typing import Dict, IO, List, Set, FrozenSet, Tuple
from moPepGen import gtf, seqvar
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from moPepGen.aa.VariantPeptideLabel import VariantPeptideInfo, \
    VariantSourceSet, LabelSourceMapping
from moPepGen.aa.VariantPeptidePool import VariantPeptidePool
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen.aa.PeptidePoolSplitter import NONCODING_SOURCE


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

class NoncanonicalPeptideSummariyTable():
    """ Summary table for noncaonincal peptides """
    KEY_N_TOTAL = 'n_total'
    def __init__(self, data:Dict[str, Dict[str,int]]=None, max_misc:int=0):
        """ constructor """
        self.data = data or {}
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

    def add_entry(self, seq:AminoAcidSeqRecord, label_map:LabelSourceMapping,
            anno:GenomicAnnotation, enzyme:str) -> None:
        """ Add peptide entry to the summary table """
        peptide_labels = VariantPeptideInfo.from_variant_peptide(
            peptide=seq, anno=anno, label_map=label_map,
            check_source=True
        )
        peptide_labels.sort()
        sources = frozenset(peptide_labels[0].sources)
        self.increment_total(sources)

        exception = 'trypsin_exception' if enzyme == 'trypsin' else None
        misc = len(seq.find_all_enzymatic_cleave_sites(enzyme, exception))
        self.increment_misc(sources, misc)

    def _increment(self, sources:FrozenSet[str], key):
        """ Increment """
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
        return ('n_total', *[self.get_key_x_misc(i) for i in range(self.max_misc)])

    def get_stringified_summary_entry(self, key:str, sep:str='\t') -> str:
        """ Get a row from the summary table for a given combination of sources. """
        rowname = '-'.join(key)
        if key not in self.data:
            entry = [rowname, '0']
            for _ in range(self.max_misc):
                entry.append('0')
        else:
            entry = [rowname, str(self.get_n_total(key))]
            for x in range(self.max_misc):
                entry.append(str(self.get_n_x_misc(key, x)))
        return sep.join(entry)

class PeptidePoolSummarizer():
    """ Summarize the variant peptide pool called by moPepGen callVariant. """
    def __init__(self, peptides:VariantPeptidePool=None,
            summary_table:NoncanonicalPeptideSummariyTable=None,
            label_map:LabelSourceMapping=None, order:Dict[str,int]=None,
            source_parser_map:Dict[str,str]=None):
        """ """
        self.peptides = peptides
        self.summary_table = summary_table or NoncanonicalPeptideSummariyTable()
        self.label_map = label_map or LabelSourceMapping()
        self.order = order or {}
        self.source_parser_map = source_parser_map or {}

    def append_order_noncoding(self):
        """ Add noncoding to the end of order """
        source = NONCODING_SOURCE
        if source not in self.order:
            self.append_order(source)

    def append_order(self, source:str):
        """ Add a group to the end of the order. """
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

    def count_peptide_source(self, anno:gtf.GenomicAnnotation, enzyme:str):
        """ Count number of peptides in each source or combinations of sources. """
        VariantSourceSet.set_levels(self.order)
        for seq in self.peptides.peptides:
            self.summary_table.add_entry(seq, self.label_map, anno, enzyme)

    def contains_exclusive_sources(self, sources:Set[str]) -> bool:
        """ Checks whether the given sources contains any mutually exclusive
        parsers. """
        for source in sources:
            if source is NONCODING_SOURCE:
                continue
            parser = self.source_parser_map[source]
            if parser not in MUTUALLY_EXCLUSIVE_PARSERS:
                continue
            # Checks if any of parsers of the rest sources are mutually
            # exclusive to the current one.
            rest_parsers = {self.source_parser_map[x] for x in sources
                if x not in {source, NONCODING_SOURCE}}
            if any(x in rest_parsers for x in MUTUALLY_EXCLUSIVE_PARSERS[parser]):
                return True
        return False

    def write_summary_table(self, handle:IO):
        """ Write summary table to output. """
        summary_keys = self.summary_table.get_keys()
        header = '\t'.join(['sources', *summary_keys])
        handle.write(header + '\n')
        sources = [it[0] for it in sorted(self.order.items(), key=lambda x:x[1])]
        for i in range(len(sources)):
            for comb in itertools.combinations(sources, i + 1):
                if self.contains_exclusive_sources(comb):
                    continue
                key = frozenset(comb)
                record = self.summary_table.get_stringified_summary_entry(key)
                handle.write(record + '\n')
