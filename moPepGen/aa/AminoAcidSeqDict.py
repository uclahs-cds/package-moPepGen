""" Model for protein or peptide sequence data
"""
from typing import Set
from Bio import SeqIO
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord


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
    
    def __setitem__(self, k:str, v:AminoAcidSeqRecord)->None:
        """ set item """
        if not isinstance(v, AminoAcidSeqRecord):
            raise ValueError('The value of a DNASeqDict must be '
            'AminoAcidSeqRecord.')
        super().__setitem__(k, v)

    def dump_fasta(self, path:str, source:str=None)->None:
        """ Dump a FASTA file into an AminoAcidSeqDict

        Args:
            path (str): Path to the FASTA file of protein sequences.
            source (str): The source of the FASTA file. Can be 'gencode' or
                'ensembl'. If None is given, it will be infered by trying
                the first 100 records. Default to None.
        """
        if not source:
            count = 0
            infered = set()
        for record in SeqIO.parse(path, 'fasta'):
            record.__class__ = AminoAcidSeqRecord
            
            if count > 100 and not source:
                source = infered.pop()
            
            if not source:
                count += 1
                infered.add(record.infer_ids(style=source))
                continue

            record.infer_ids(source)
            if record.transcript_id in self.keys():
                raise ValueError(
                    'Duplicated seqnames found in FASTA file: ' + path
                )
            self[record.transcript_id] = record

    def create_unique_peptide_pool(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500.)->Set[str]:
        """ Create a unique piptide pool. All peptides with I (isoleucine)
        replaced with L (leucine) are also added to the pool. This is used to
        filter variant peptides. For proteins that starts with X, the part of
        the sequence before the first cleave site is removed. Protein sequences
        that starts with X usually means that the 5' CDS is incomplete.

        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            start (int): Index to start searching.
            min_mw (float): Minimal molecular weight of the peptides to report.
                Defaults to 500.
        
        Returns:
            A set of unique peptides as string.
        """
        pool = set()
        protein: AminoAcidSeqRecord
        for protein in self.values():
            if protein.seq.startswith('X'):
                try:
                    cleave_site = protein.find_first_enzymatic_cleave_sites(
                        rule=rule, exception=exception)
                except StopIteration:
                    continue
                protein = protein[cleave_site:]
            peptides = protein.enzymatic_cleave(rule=rule, exception=exception,
                miscleavage=miscleavage, min_mw=min_mw)
            for peptide in peptides:
                pool.add(str(peptide.seq))
                # Convert all I with L and add to the carnonical peptide pool.
                pool.add(str(peptide.seq).replace('I', 'L'))
        return pool


if __name__ == '__main__':
    proteome = AminoAcidSeqDict()
    file_dir = 'test/files/downsampled_set'
    proteome.dump_fasta(f'{file_dir}/gencode_v34_translations_chr22.fasta')
    proteome.create_unique_peptide_pool('trypsin')