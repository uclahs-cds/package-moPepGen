""" module for peptide peptide labels """
from __future__ import annotations
from typing import Dict, Iterable, List, TYPE_CHECKING
from moPepGen import err, seqvar, circ,SPLIT_DATABASE_KEY_SEPARATER
from . import VariantPeptideIdentifier as pi

if TYPE_CHECKING:
    from .AminoAcidSeqRecord import AminoAcidSeqRecord
    from moPepGen.gtf import GenomicAnnotation

NONCODING_SOURCE = 'Noncoding'

class VariantSourceSet(set):
    """ Variant source set. This is a class of ordered set.

    Example:
        >>> VariantSourceSet.set_levels({'A':0, 'B':1, 'C':2})
        >>> print(VariantSourceSet(['A']) < VariantSourceSet(['B']))
        True

        >>> print(VariantSourceSet(['A', 'B']) > VariantSourceSet(['B']))
        True
    """
    levels_map = {}
    levels = []

    def __init__(self, *args):
        """ constructor """
        for it in list(*args):
            self.validate(it)
        super().__init__(*args)

    def validate(self, it):
        """ Validate element """
        if it not in self.levels_map:
            raise ValueError(f'No defined level found for {it}')

    def add(self, element):
        """ override the add method """
        self.validate(element)
        super().add(element)

    @classmethod
    def set_levels(cls, levels:Dict[str,int]):
        """ set levels. This is to make sure that all instances share the same
        order. """
        cls.levels_map = levels
        cls.levels = [y[0] for y in sorted(cls.levels_map.items(), key=lambda x:x[1])]

    @classmethod
    def reset_levels(cls):
        """ reset levels """
        cls.levels_map = {}
        cls.levels = []

    def __str__(self) -> str:
        """ str """
        sorted_list = [x for x in self.levels if x in self]
        return SPLIT_DATABASE_KEY_SEPARATER.join(sorted_list)

    def __gt__(self, other:VariantSourceSet) -> bool:
        """ greater than """
        if self == other:
            return False
        if len(self) > len(other):
            return True
        if len(self) < len(other):
            return False
        this = self.to_int()
        that = other.to_int()
        for i, j in zip(this, that):
            if i > j:
                return True
            if i < j:
                return False
        return False

    def __ge__(self, other:VariantSourceSet) -> bool:
        """ greater or equal to """
        return self == other or self > other

    def __lt__(self, other:VariantSourceSet) -> bool:
        """ less then """
        return not self >= other

    def __le__(self, other:VariantSourceSet) -> bool:
        """ less or equal to """
        return not self > other

    def to_int(self, sort=True) -> Iterable[int]:
        """ to int """
        source_int = {self.levels_map[x] for x in self}
        if sort:
            source_int = list(source_int)
            source_int.sort()
        return source_int

def is_circ_rna(_id:str) -> bool:
    """ Check if the id is a circRNA """
    return _id.startswith('CIRC-') or _id.startswith('CI-')

class VariantPeptideInfo():
    """ Variant peptide label. This is a helper class in order to sort peptide
    labels easily.

    Args:
        transcript_id (str): If the variant is a circRNA, this is the circRNA
            ID.
    """
    def __init__(self, orignial_label:str, gene_ids:List[str],
            variant_labels:Dict[str,List[str]],
            variant_index:int, sources:VariantSourceSet=None):
        """ Constructor """
        self.original_label = orignial_label
        self.gene_ids = gene_ids
        self.variant_labels = variant_labels
        self.variant_index = variant_index
        self.sources = sources or VariantSourceSet()


    @staticmethod
    def from_variant_peptide_minimal(peptide:AminoAcidSeqRecord
            ) -> List[VariantPeptideInfo]:
        """ Create list of VariantPeptideInfo with minimal information. """
        info_list:List[VariantPeptideInfo] = []
        variant_ids = pi.parse_variant_peptide_id(peptide.description)
        for variant_id in variant_ids:
            if isinstance(variant_id, pi.CircRNAVariantPeptideIdentifier):
                gene_ids = None
                var_ids = {}

            elif isinstance(variant_id, pi.FusionVariantPeptideIdentifier):
                first_tx_id = variant_id.first_tx_id
                second_tx_id = variant_id.second_tx_id
                gene_ids = None
                var_ids = {
                    first_tx_id: variant_id.first_variants + [variant_id.fusion_id],
                    second_tx_id: variant_id.second_variants
                }

            elif isinstance(variant_id, pi.BaseVariantPeptideIdentifier):
                gene_ids = None
                var_ids = {}

            info = VariantPeptideInfo(str(variant_id), gene_ids, var_ids, variant_id.index)

            info_list.append(info)
        return info_list

    @staticmethod
    def from_variant_peptide(peptide:AminoAcidSeqRecord,
            anno:GenomicAnnotation, label_map:LabelSourceMapping=None,
            check_source:bool=True
            ) -> List[VariantPeptideInfo]:
        """ Parse from a variant peptide record """
        info_list = []
        variant_ids = pi.parse_variant_peptide_id(peptide.description)
        for variant_id in variant_ids:
            if isinstance(variant_id, pi.CircRNAVariantPeptideIdentifier):
                circ_rna_id = variant_id.circ_rna_id
                tx_id = circ_rna_id.split('-', 2)[1]
                gene_ids = [anno.transcripts[tx_id].transcript.gene_id]
                var_ids = {gene_ids[0]: [circ_rna_id, *variant_id.variant_ids]}

            elif isinstance(variant_id, pi.FusionVariantPeptideIdentifier):
                first_tx_id = variant_id.first_tx_id
                first_gene_id = anno.transcripts[first_tx_id].transcript.gene_id
                second_tx_id = variant_id.second_tx_id
                second_gene_id = anno.transcripts[second_tx_id].transcript.gene_id
                gene_ids = [first_gene_id, second_gene_id]
                var_ids = {
                    first_gene_id: variant_id.first_variants + [variant_id.fusion_id],
                    second_gene_id: variant_id.second_variants
                }

            elif isinstance(variant_id, pi.BaseVariantPeptideIdentifier):
                gene_id = anno.transcripts[variant_id.transcript_id].transcript.gene_id
                gene_ids = [gene_id]
                var_ids = {gene_id: variant_id.variant_ids}

            info = VariantPeptideInfo(str(variant_id), gene_ids, var_ids, variant_id.index)

            if check_source:
                has_noncoding = False
                for gene_id in gene_ids:
                    gene_model = anno.genes[gene_id]
                    for tx_id in gene_model.transcripts:
                        if not anno.transcripts[tx_id].is_protein_coding:
                            info.sources.add(NONCODING_SOURCE)
                            has_noncoding = True
                            break
                    if has_noncoding:
                        break

                for gene_id, _ids in var_ids.items():
                    for var_id in _ids:
                        source = label_map.get_source(gene_id, var_id)
                        info.sources.add(source)

            info_list.append(info)
        return info_list

    def is_fusion(self) -> bool:
        """ Check if this is a fusion """
        _id = pi.parse_variant_peptide_id(self.original_label)[0]
        return isinstance(_id, pi.FusionVariantPeptideIdentifier)

    def is_circ_rna(self) -> bool:
        """ Check if this is a circRNA """
        _id = pi.parse_variant_peptide_id(self.original_label)[0]
        return isinstance(_id, pi.CircRNAVariantPeptideIdentifier)

    def is_splice_altering(self) -> bool:
        """ Check if the variant paptide label is alternative splicing """
        _id = pi.parse_variant_peptide_id(self.original_label)[0]
        return isinstance(_id, pi.BaseVariantPeptideIdentifier) and \
            _id.is_alternative_splicing()

    @staticmethod
    def is_noncoding(tx_id:str, inclusion:List[str], exclusion:List[str]) -> bool:
        """ Check if a transcript is noncoding """
        if inclusion:
            return tx_id in inclusion
        if exclusion:
            return tx_id not in exclusion
        raise ValueError('At least one of inclusion or exclusion must given')


    def all_noncoding(self, anno:GenomicAnnotation) -> bool:
        """ Check if all transcripts are noncoding """
        tx_ids = self.get_transcript_ids()
        return all(not anno.transcripts[tx_id].is_protein_coding for tx_id in tx_ids)

    def any_noncoding(self, anno:GenomicAnnotation) -> bool:
        """ Check if all transcripts are noncoding """
        tx_ids = self.get_transcript_ids()
        return any(not anno.transcripts[tx_id].is_protein_coding for tx_id in tx_ids)

    def all_coding(self, anno:GenomicAnnotation) -> bool:
        """ Check if all transcripts are coding """
        tx_ids = self.get_transcript_ids()
        return all(anno.transcripts[tx_id].is_protein_coding for tx_id in tx_ids)

    def any_coding(self, anno:GenomicAnnotation) -> bool:
        """ Check if all transcripts are coding """
        tx_ids = self.get_transcript_ids()
        return any(anno.transcripts[tx_id].is_protein_coding for tx_id in tx_ids)

    def get_transcript_ids(self) -> List[str]:
        """ get transcript IDs """
        variant_id = pi.parse_variant_peptide_id(self.original_label)[0]
        if isinstance(variant_id, pi.CircRNAVariantPeptideIdentifier):
            return [variant_id.circ_rna_id.split('-', 2)[1]]
        if isinstance(variant_id, pi.FusionVariantPeptideIdentifier):
            return [variant_id.first_tx_id, variant_id.second_tx_id]
        if isinstance(variant_id, pi.BaseVariantPeptideIdentifier):
            return [variant_id.transcript_id]
        raise ValueError('Variant ID unrecognized')

    def __str__(self) -> str:
        """ str """
        return self.original_label

    def __eq__(self, other:VariantPeptideInfo):
        """ equal to """
        if len(self.gene_ids) != len(other.gene_ids):
            return False
        if any(x != y for x,y in zip(self.gene_ids, other.gene_ids)):
            return False
        return set(self.variant_labels) == set(other.variant_labels) \
            and self.variant_index == other.variant_index \
            and self.sources == other.sources

    def __ne__(self, other:VariantPeptideInfo):
        """ not equal to """
        return not self == other

    def __gt__(self, other:VariantPeptideInfo):
        """" greater than """
        return self.sources > other.sources

    def __ge__(self, other:VariantPeptideInfo):
        """ greater than or equal to """
        return self > other or self == other

    def __lt__(self, other:VariantPeptideInfo):
        """ less than """
        return not self >=  other

    def __le__(self, other:VariantPeptideInfo):
        """ less than or equal to """
        return not self > other

class LabelSourceMapping():
    """ Helper class to handle label source mapping """
    def __init__(self, data:Dict[str,Dict[str,str]]=None):
        """ construnctor"""
        self.data = data or {}

    def add_variant(self, variant:seqvar.VariantRecord, source:str):
        """ add variant """
        gene_id = variant.location.seqname
        if gene_id not in self.data:
            self.data[gene_id] = {}
        if variant.id not in self.data[gene_id]:
            self.data[gene_id][variant.id] = source

    def add_circ_rna(self, record:circ.CircRNAModel, source:str):
        """ add circRNA record """
        gene_id = record.gene_id
        if gene_id not in self.data:
            self.data[gene_id] = {}
        if record.id not in self.data[gene_id]:
            self.data[gene_id][record.id] = source

    def get_source(self, gene_id:str, var_id:str) -> bool:
        """ Get source """
        try:
            return self.data[gene_id][var_id]
        except KeyError as e:
            raise err.VariantSourceNotFoundError(gene_id, var_id) from e
