""" Module for Amino Acid Record """
from __future__ import annotations
from typing import Iterable, List, Set
import re
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from moPepGen.aa.expasy_rules import EXPASY_RULES


class AminoAcidSeqRecord(SeqRecord):
    #pylint: disable=W0223
    """ A AminoAcidSeqRecord holds a protein or peptide sequence
    """
    def __init__(self, seq:SeqRecord, _id:str="<unknown id>",
            name:str="<unknown name>", description:str="<unknown description>",
            gene_id:str=None, transcript_id:str=None, protein_id:str=None,
            gene_name:str=None, **kwargs):
        self.id = None
        self.name = None
        super().__init__(seq, id=_id, name=name, description=description,
            **kwargs)
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.protein_id = protein_id
        self.gene_name = gene_name

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
        """ Equal to. It is implemented in this way, because the get_equivalent
        calls this __eq__ instead of the other for reason that I don't
        understand. """
        result = (self.seq == other.seq)
        if result and hasattr(other, 'match'):
            other.match = self
        return result

    def __ne__(self, other:AminoAcidSeqRecord) -> bool:
        """ not equal to """
        return not self == other

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

        self.gene_id = gene_id.split('.')[0]
        self.protein_id = protein_id.split('.')[0]
        self.transcript_id = transcript_id.split('.')[0]

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
            miscleavage:int=2, min_mw:float=500.0, min_length:int=7,
            max_length:int=25, cds_start_nf:bool=False
            )->Set[AminoAcidSeqRecord]:
        """ Performs enzymatic cleave """
        peptides = []
        sites = [0]
        sites += self.find_all_enzymatic_cleave_sites(
            rule=rule, exception=exception)
        sites.append(len(self))
        start = 0

        def update_peptides(peptide):
            mol_wt = SeqUtils.molecular_weight(peptide.seq, 'protein')
            weight_flag = mol_wt > min_mw
            length_flag = len(peptide.seq) >= min_length \
                and len(peptide.seq) <= max_length
            if weight_flag and length_flag:
                peptides.append(peptide)

        while start < len(sites) - 1:
            end = start + 1
            while end - start - 1 <= miscleavage and end < len(sites):
                peptide = self[sites[start]:sites[end]]
                if start == 0 and not cds_start_nf and peptide.seq.startswith('M'):
                    update_peptides(peptide[1:])
                update_peptides(peptide)
                end += 1
            start += 1
        return peptides
