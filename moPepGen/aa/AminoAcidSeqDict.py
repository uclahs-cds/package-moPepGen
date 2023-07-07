""" Model for protein or peptide sequence data
"""
from __future__ import annotations
from typing import Set, TYPE_CHECKING
from Bio import SeqIO
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from moPepGen.version import MetaVersion


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation

class AminoAcidSeqDict(dict):
    """ AminoAcidSeqDict """
    def __init__(self, *args, **kwargs):
        """ Constructor """
        for val in kwargs.values():
            if not isinstance(val, AminoAcidSeqRecord):
                raise ValueError(
                    'The value of a DNASeqDict must be AminoAcidSeqRecord.'
                )
        super().__init__(*args, **kwargs)
        self.version = MetaVersion()

    def __setitem__(self, k:str, v:AminoAcidSeqRecord)->None:
        """ set item """
        if not isinstance(v, AminoAcidSeqRecord):
            raise ValueError('The value of a AASeqDict must be '
            'AminoAcidSeqRecord.')
        super().__setitem__(k, v)

    def dump_fasta(self, path:str, source:str=None)->None:
        """ Dump a FASTA file into an AminoAcidSeqDict

        Args:
            path (str): Path to the FASTA file of protein sequences.
            source (str): The source of the FASTA file. Can be 'GENCODE' or
                'ENSEMBL'. If None is given, it will be infered by trying
                the first 100 records. Default to None.
        """
        count = 0
        if not source:
            inferred = set()
        record:AminoAcidSeqRecord
        for record in SeqIO.parse(path, 'fasta'):
            record.__class__ = AminoAcidSeqRecord
            if count > 100 and not source:
                source = inferred.pop()

            if not source:
                count += 1
                inferred.add(record.infer_ids(style=source))
            else:
                record.infer_ids(source)
            if record.transcript_id in self.keys():
                raise ValueError(
                    'Duplicated seqnames found in FASTA file: ' + path
                )
            self[record.transcript_id] = record

    def create_unique_peptide_pool(self, anno:GenomicAnnotation,
            rule:str, exception:str=None, miscleavage:int=2, min_mw:float=500.,
            min_length:int=7, max_length:int=25)->Set[str]:
        """ Create a unique piptide pool.

        All peptides with I (isoleucine)
        replaced with L (leucine) are also added to the pool. This is used to
        filter variant peptides. For proteins that starts with X, the part of
        the sequence before the first cleave site is removed. Protein sequences
        that starts with X usually means that the 5' CDS is incomplete.

        For a given gene, if the cds start is known (i.e. the gene model does not
        have the cds_start_NF tag from the GTF file), the peptide with and
        without the leading methionine will both be included. This is because
        sometimes the start codon methionine is cleaved spontaneously inside the
        cell.

        Args:
            anno (GenomicAnnotation): Genomic annotation parsed from GTF.
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            start (int): Index to start searching.
            min_mw (float): Minimal molecular weight of the peptides to report.
                Defaults to 500.
            min_length (int): Minimal length of the peptides to report, inclusive.
                Defaults to 7.
            max_length (int): Maximum length of the peptides to report, inclusive.
                Defaults to 25.

        Returns:
            A set of unique peptides as string.
        """
        pool = set()
        protein: AminoAcidSeqRecord
        it = iter(self.values())
        protein = next(it, None)
        while protein:
            tx_id = protein.transcript_id
            if tx_id in anno.transcripts:
                cds_start_nf = anno.transcripts[tx_id].is_cds_start_nf()
            else:
                cds_start_nf = False
            if protein.seq.startswith('X'):
                protein.seq = protein.seq.lstrip('X')

            stop_site = protein.seq.find('*')
            if stop_site > -1:
                protein = protein[:stop_site]

            try:
                peptides = protein.enzymatic_cleave(
                    rule=rule,
                    exception=exception,
                    miscleavage=miscleavage,
                    min_mw=min_mw,
                    min_length=min_length,
                    max_length=max_length,
                    cds_start_nf=cds_start_nf
                )
            except ValueError as e:
                msg = "'X' is not a valid unambiguous letter for protein"
                if e.args[0] == msg:
                    protein.seq = protein.seq.split('X')[0]
                    continue
                raise e
            for peptide in peptides:
                pool.add(str(peptide.seq))
                # Convert all I with L and add to the canonical peptide pool.
                pool.add(str(peptide.seq).replace('I', 'L'))

            protein = next(it, None)

        return pool
