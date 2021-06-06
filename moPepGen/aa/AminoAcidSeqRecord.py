""" Module for Amino Acid Record """
from __future__ import annotations
from typing import Iterable, List, Set
import re
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from moPepGen.aa.expasy_rules import EXPASY_RULES


_NO_SEQRECORD_COMPARISON = "SeqRecord comparison is deliberately not" +\
    " implemented. Explicitly compare the attributes of interest."
class AminoAcidSeqRecord(SeqRecord):
    """ A AminoAcidSeqRecord holds a protein or peptide sequence
    """
    def __init__(self, seq:SeqRecord, _id:str="<unknown id>",
            name:str="<unknown name>", description:str="<unknown description>",
            gene_id:str=None, transcript_id:str=None, protein_id:str=None,
            **kwargs):
        super().__init__(seq, id=_id, name=name, description=description,
            **kwargs)
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.protein_id = protein_id

    def __getitem__(self, index) -> AminoAcidSeqRecord:
        """"""
        new_seq = self.seq.__getitem__(index)
        new_one = self.__class__(seq=new_seq, _id=self.id, name=self.name,
            description=self.description, gene_id=self.gene_id,
            transcript_id=self.transcript_id, protein_id=self.protein_id)
        return new_one

    def __add__(self, other:AminoAcidSeqRecord) -> AminoAcidSeqRecord:
        """"""
        new_one = super().__add__(other)
        new_one.__class__ = self.__class__
        new_one.gene_id = self.gene_id
        new_one.transcript_id = self.transcript_id
        new_one.protein_id = self.protein_id
        return new_one

    def __hash__(self):
        """ hash """
        return hash(str(self.seq))

    def __eq__(self, other:AminoAcidSeqRecord):
        """ equal to """
        return self.seq == other.seq

    def __ne__(self, other:AminoAcidSeqRecord) -> bool:
        """ not equal to """
        return not self == other

    def __gt__(self, other:AminoAcidSeqRecord):
        """ greater than """
        return NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __ge__(self, other:AminoAcidSeqRecord):
        """ greater or equal to """
        return NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __lt__(self, other:AminoAcidSeqRecord):
        """ less than """
        return NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __le__(self, other:AminoAcidSeqRecord):
        """ less or equal to """
        return NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __le___(self, other:AminoAcidSeqRecord):
        """ Due to an typo in biopython """
        return NotImplementedError(_NO_SEQRECORD_COMPARISON)

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
                self.infer_ids_gencode()
                return 'gencode'
            except ValueError:
                pass

            try:
                self.infer_ids_ensembl()
                return 'ensembl'
            except ValueError as e:
                raise ValueError(
                    'Failed to infer gene ID, transcript ID, and protein ID '
                    r'using both ENSEMBL and GENCODE\'s format.'
                ) from e
        elif style == 'gencode':
            self.infer_ids_gencode()
        elif style == 'ensembl':
            self.infer_ids_ensembl()
        else:
            raise ValueError(f'style {style} is not supported')
        return style


    def infer_ids_gencode(self):
        """ Infers the gene, transcript, and protein ID from description base
        on GENCODE's format
        """
        delimiter = '|'
        if delimiter not in self.description:
            raise ValueError(r'The description does not have any \'' +
            delimiter + r'\'.')
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
            raise ValueError(r'Couldn\'t find gene ID')
        if transcript_id is None:
            raise ValueError(r'Transcript ID could\'t be inferred.')
        if protein_id is None:
            raise ValueError(r'Protein ID could\'t be inferred.')

        self.id = protein_id
        self.name = protein_id
        self.gene_id = gene_id
        self.protein_id = protein_id
        self.transcript_id = transcript_id

    def infer_ids_ensembl(self):
        """ Infers the gene, transcript, and protein ID from description base
        on ENSEMBL's format
        """
        if not re.search(r'^ENSP\d+\.\d+ pep ', self.description):
            raise ValueError(r'Description doesn\'t seem like ESEMBL\'s style')

        protein_id = self.id

        match = re.search(r'gene:(ENSG\d+\.\d+) ', self.description)
        if not match:
            raise ValueError(r'Couldn\'t find gene ID')
        gene_id = match.group(1)

        match = re.search(r'transcript:(ENST\d+\.\d+) ', self.description)
        if not match:
            raise ValueError(r'Couldn\'t find transcript ID')
        transcript_id = match.group(1)

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

        for it in re.finditer(rule, str(self.seq)):
            if it not in exceptions:
                yield it.end()

    def find_first_enzymatic_cleave_site(self, rule:str, exception:str=None,
            start:int=0) -> int:
        """ Find the first enzymatic cleave site """
        it = self[start:]\
            .iter_enzymatic_cleave_sites(rule=rule, exception=exception)
        try:
            return next(it) + start
        except StopIteration:
            return -1

    def find_all_enzymatic_cleave_sites(self, rule:str, exception:str=None,
            ) -> List[int]:
        """ Find all enzymatic cleave sites. """
        return list(self.iter_enzymatic_cleave_sites(rule=rule,
            exception=exception))

    def enzymatic_cleave(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500.0)->Set[AminoAcidSeqRecord]:
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
