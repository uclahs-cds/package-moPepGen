""""""
from __future__ import annotations
from typing import Iterable
import re
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from moPepGen.aa.expasy_rules import EXPASY_RULES


class AminoAcidSeqRecord(SeqRecord):
    """ A AminoAcidSeqRecord holds a protein or peptide sequence
    """
    def __init__(self, seq:SeqRecord, id:str="<unknown id>",
            name:str="<unknown name>", description:str="<unknown description>",
            gene_id:str=None, transcript_id:str=None, protein_id:str=None,
            **kwargs):
        super().__init__(seq, id=id, name=name, description=description,
            **kwargs)
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.protein_id = protein_id
    
    def __hash__(self):
        """ hash """
        return hash(str(self.seq))
    
    def __eq__(self, other:AminoAcidSeqRecord):
        """ equal to """
        return self.seq == other.seq
    
    def infer_ids(self, style:str=None) -> str:
        """
        Args:
            source (str): The style of the fasta file header. Either 'gencode'
                or 'ensembl'. If None, it will try both. Defaults to None.
        
        Returns:
            The style ifered.
        """
        if not style:
            try:
                self._infer_ids_gencode()
                return 'gencode'
            except ValueError:
                pass
            
            try:
                self._infer_ids_ensembl()
                return 'ensembl'
            except ValueError:
                raise ValueError(
                    'Failed to infer gene ID, transcript ID, and protein ID '
                    'using both ENSEMBL and GENCODE\'s format.'
                )
        elif style == 'gencode':
            self._infer_ids_gencode()
        elif style == 'ensembl':
            self._infer_ids_ensembl()
        else:
            raise ValueError(f'style {style} is not supported')
        return style


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

    def iter_enzymatic_cleave_sites(self, rule:str, exception:str=None
            ) -> Iterable[int]:
        """ Create a generator of the cleave sites """
        rule = EXPASY_RULES[rule]
        exception = EXPASY_RULES.get(exception, exception)
        exceptions = [] if exception is None else \
            [x.end() for x in re.finditer(exception, str(self.seq))]
        
        for x in re.finditer(rule, str(self.seq)):
            if x not in exceptions:
                yield x.end()
        return
    
    def find_first_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            start:int=0) -> int:
        """ Find the first enzymatic cleave site """
        iter = self[start:]\
            .iter_enzymatic_cleave_sites(rule=rule, exception=exception)
        return next(iter) + start
    
    def find_all_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            ) -> List[int]:
        """"""
        return [i for i in self.iter_enzymatic_cleave_sites(rule=rule,
            exception=exception)]

    def enzymatic_cleave(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500.0)->set[AminoAcidSeqRecord]:
        """ Performs enzymatic cleave """
        peptides = []
        sites = [0]
        sites += self.find_all_enzymatic_cleave_sites(
            rule=rule, exception=exception)
        sites.append(len(self))
        start = 0
        while start < len(sites) - 1:
            end = start + 1
            while end - start - 1 <= miscleavage and end < len(sites):
                peptide = self[sites[start]:sites[end]]
                if SeqUtils.molecular_weight(peptide.seq, 'protein') > min_mw:
                    peptides.append(peptide)
                end += 1
            start += 1
        return peptides
        