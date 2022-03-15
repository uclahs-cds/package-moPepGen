""" module for peptide pool summarizer """
from __future__ import annotations
import itertools
from typing import Dict, IO
from Bio import SeqIO
from moPepGen import gtf, seqvar, aa
from moPepGen.aa.VariantPeptideLabel import VariantPeptideInfo, \
    VariantSourceSet, LabelSourceMapping
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen.aa.PeptidePoolSplitter import NONCODING_SOURCE

class PeptidePoolSummarizer():
    """ Summarize the variant peptide pool called by moPepGen callVariant. """
    def __init__(self, summary_table:Dict[VariantSourceSet,int]=None,
            label_map:LabelSourceMapping=None, order:Dict[str,int]=None):
        """ """
        self.summary_table = summary_table or {}
        self.label_map = label_map or LabelSourceMapping()
        self.order = order or {}

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

    def update_label_map(self, handle:IO):
        """ Read GVF files and add to the label source map. """
        metadata = GVFMetadata.parse(handle)
        self.append_order(metadata.source)
        for gene_id, _, label in seqvar.io.parse_label(handle):
            self.label_map.add_record(gene_id, label, metadata.source)

    def count_peptide_source(self, handle:IO, anno:gtf.GenomicAnnotation):
        """ Count number of peptides in each source or combinations of sources. """
        VariantSourceSet.set_levels(self.order)
        for seq in SeqIO.parse(handle, 'fasta'):
            seq.__class__ = aa.AminoAcidSeqRecord
            seq.id = seq.description
            seq.name = seq.description
            peptide_labels = VariantPeptideInfo.from_variant_peptide(
                peptide=seq, anno=anno, label_map=self.label_map,
                check_source=True
            )
            peptide_labels.sort()
            sources = frozenset(peptide_labels[0].sources)
            if sources not in self.summary_table:
                self.summary_table[sources] = 1
            else:
                self.summary_table[sources] += 1

    def write_summary_table(self, handle:IO):
        """ Write summary table to output. """
        handle.write("sources\tn_peptides\n")
        sources = [it[0] for it in sorted(self.order.items(), key=lambda x:x[1])]
        for i in range(len(sources)):
            for comb in itertools.combinations(sources, i + 1):
                key = frozenset(comb)
                if key in self.summary_table:
                    count = self.summary_table[key]
                else:
                    count = 0
                handle.write(f"{','.join(comb)}\t{count}\n")
