""" Module for variant peptide pool (unique) """
from __future__ import annotations
from typing import Set, IO, Dict, List, Tuple
from pathlib import Path
from Bio import SeqUtils, SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from moPepGen import get_equivalent, VARIANT_PEPTIDE_SOURCE_DELIMITER
from .VariantPeptideLabel import VariantPeptideInfo


class VariantPeptidePool():
    """ Varaint Peptide Pool """
    def __init__(self, peptides:Set[AminoAcidSeqRecord]=None):
        """ Constructor """
        self.peptides = peptides or set()
        self.peptide_delimeter = VARIANT_PEPTIDE_SOURCE_DELIMITER

    def add_peptide(self, peptide:AminoAcidSeqRecord,
            canonical_peptides:Set[str], min_mw:int=500, min_length:int=7,
            max_length:int=25, skip_checking:bool=False):
        """ Add a peptide to the pool if it does not already exist. Otherwise,
        the label is appended to the existing same peptide.

        Args:
            peptide (AminoAcidSeqRecord): The amino acid sequence to be added.
            min_mw (int): Minimal molecular weight.
            min_length (int): Minimal peptide sequence length.
            max_length (int): Maximal peptide sequence length.
            canonical_peptides (Set[str]): Canonical peptides.
        """
        if not skip_checking:
            if SeqUtils.molecular_weight(peptide.seq, 'protein') < min_mw:
                return
            if len(peptide.seq) < min_length or len(peptide.seq) > max_length:
                return
            if str(peptide.seq) in canonical_peptides:
                return
        same_peptide = get_equivalent(self.peptides, peptide)
        if same_peptide:
            same_peptide:Seq
            new_label = peptide.description
            same_peptide.description += (self.peptide_delimeter + new_label)
            same_peptide.id = same_peptide.description
            same_peptide.name = same_peptide.description
        else:
            self.peptides.add(peptide)

    def write(self, path:Path):
        """ Write the variant peptide pool to FASTA file.

        Args:
            path (Path): Output FASTA file path.
        """
        with open(path, 'w') as handle:
            record2title = lambda x: x.description
            writer = FastaIO.FastaWriter(handle, record2title=record2title)
            for record in self.peptides:
                writer.write_record(record)

    @classmethod
    def load(cls, handle:IO) -> VariantPeptidePool:
        """ Load variant peptide pool from file """
        pool = cls()
        for seq in SeqIO.parse(handle, 'fasta'):
            seq.__class__ = AminoAcidSeqRecord
            seq.id = seq.description
            seq.name = seq.description
            pool.peptides.add(seq)
        return pool

    def filter(self, exprs:Dict[str,int]=None, cutoff:float=None,
            coding_transcripts:List[str]=None, keep_all_noncoding:bool=False,
            keep_all_coding:bool=False, enzyme:str='trypsin',
            miscleavage_range:Tuple[int,int]=(None, None)
            ) -> VariantPeptidePool:
        """ Filter variant peptides according to gene expression. """
        label_delimiter = VARIANT_PEPTIDE_SOURCE_DELIMITER
        filtered_pool = VariantPeptidePool()
        for peptide in self.peptides:
            if any(x is not None for x in miscleavage_range):
                exception = 'trypsin_exception' if enzyme == 'trypsin' else None
                misc = peptide.find_all_enzymatic_cleave_sites(enzyme, exception)
                if miscleavage_range[0] is not None and len(misc) < miscleavage_range[0]:
                    continue
                if miscleavage_range[1] is not None and len(misc) > miscleavage_range[1]:
                    continue
            peptide_entries = VariantPeptideInfo.from_variant_peptide_minimal(peptide)
            keep = []
            for entry in peptide_entries:
                all_noncoding = not any(x in coding_transcripts
                    for x in entry.get_transcript_ids())
                all_coding = all(x in coding_transcripts
                    for x in entry.get_transcript_ids())

                if keep_all_noncoding and all_noncoding:
                    should_keep = True

                elif keep_all_coding and all_coding:
                    should_keep = True
                else:
                    if exprs is not None:
                        tx_ids = entry.get_transcript_ids()
                        should_keep = entry.is_fusion() or entry.is_circ_rna() or\
                            entry.is_splice_altering() or \
                            all(exprs[tx] > cutoff for tx in tx_ids)
                    else:
                        should_keep = True
                if should_keep:
                    keep.append(entry)
            if keep:
                label = label_delimiter.join([str(x) for x in keep])
                peptide.description = label
                filtered_pool.peptides.add(peptide)
        return filtered_pool
