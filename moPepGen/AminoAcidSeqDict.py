""" Model for protein or peptide sequence data
"""
from __future__ import annotations
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class AminoAcidSeqRecord(SeqRecord):
    """ A AminoAcidSeqRecord holds a protein or peptide sequence
    """
    def __init__(self, seq:SeqRecord, id:str, name:str, gene_id:str,
            transcript_id:str, protein_id:str, **kwargs):
            super().__init__(
                seq, id=id, name=name, **kwargs
            )
            self.gene_id = gene_id
            self.transcript_id = transcript_id
            self.protein_id = protein_id
    
    def infer_ids(self):
        """
        """
        try:
            self._infer_ids_from_description_gencode()
            return
        except ValueError:
            pass
        
        try:
            self._infer_ids_from_description_ensembl()
        except ValueError:
            raise ValueError(
                'Failed to infer gene ID, transcript ID, and protein ID using '
                'both ENSEMBL and GENCODE\'s format.'
            )


    def _infer_ids_gencode(self):
        """ Infers the gene, transcript, and protein ID from description base
        on GENCODE's format
        """
        delimiter = '|'
        if delimiter not in self.description:
            raise ValueError('The description does not have any \'' + delimiter
            + '\'.')
        ids = self.description.split(delimiter)
        gene_id, protein_id, transcript_id = None, None, None
        for _id in ids:
            if _id.startswith('ENSP'):
                protein_id = _id
            elif _id.startswith('ENSG'):
                gene_id = _id
            elif _id.startswith('ENST'):
                transcript_id = _id
            else:
                continue
        if gene_id is None:
            raise ValueError('Couldn\'t find gene ID')
        if transcript_id is None:
            raise ValueError('Transcript ID could\'t be inferred.')
        if protein_id is None:
            raise ValueError('Protein ID could\'t be inferred.')
        
        self.gene_id = gene_id
        self.protein_id = protein_id
        self.transcript_id = transcript_id
    
    def _infer_ids_ensembl(self):
        """ Infers the gene, transcript, and protein ID from description base
        on ENSEMBL's format
        """
        if not re.search('^ENSP\d+\.\d+ pep ', self.description):
            raise ValueError('Description doesn\'t seem like ESEMBL\'s style')
        
        protein_id = self.id

        match = re.search('gene:(ENSG\d+\.\d+) ', self.description)
        if not match:
            raise ValueError('Couldn\'t find gene ID')
        gene_id = match.group(0)

        match = re.search('transcript:(ENST\d+\.\d+) ')
        if not match:
            raise ValueError('Couldn\'t find transcript ID')
        transcript_id = match.group(0)
        
        self.gene_id = gene_id
        self.protein_id = protein_id
        self.transcript_id = transcript_id


    def digest(self, enzyme)->set[AminoAcidSeqRecord]:
        pass

class AminoAcidSeqDict(dict):
    """"""
    def __init__(self, *args, **kwargs):
        for val in kwargs.values():
            if not isinstance(val, AminoAcidSeqRecord):
                raise ValueError(
                    'The value of a DNASeqDict must be AminoAcidSeqRecord.'
                )
        super().__init__(*args, **kwargs)
    
    def __setitem__(self, k:str, v:AminoAcidSeqRecord)->None:
        if not isinstance(v, AminoAcidSeqRecord):
            raise ValueError('The value of a DNASeqDict must be '
            'AminoAcidSeqRecord.')
        super().__setitem__(k, v)

    def dump_fasta(self, path:str)->None:
        """ Dump a FASTA file into an AminoAcidSeqDict

        Args:
            path (str): Path to the FASTA file of protein sequences.
        """
        for record in SeqIO.parse(path, 'fasta'):
            record.__class__ = AminoAcidSeqRecord
            record._infer_ensembl_ids_from_description()
            if record.id in self.keys():
                raise ValueError(
                    'Duplicated seqnames found in FASTA file: ' + path
                )
            self[record.id] = record

    def call_unique_peptides(self)->set[AminoAcidSeqRecord]:
        pass
    