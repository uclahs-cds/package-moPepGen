""" Module for VariantPeptideIdentifier """
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Dict, List, TYPE_CHECKING
from moPepGen import VARIANT_PEPTIDE_SOURCE_DELIMITER

if TYPE_CHECKING:
    from moPepGen.seqvar import VariantRecord

def create_variant_peptide_id(transcript_id:str, variants:List[VariantRecord],
        orf_id:str=None, index:int=None) -> str:
    """ Create variant peptide ID """
    variant_id_map:Dict[str,List[VariantRecord]] = {}
    is_fusion = False
    is_circ_rna = False
    for variant in variants:
        if variant.is_fusion():
            is_fusion = True
            fusion_variant = variant
        elif variant.is_circ_rna():
            is_circ_rna = True
            circ_rna_id = variant.id
        else:
            seqname = variant.location.seqname
            if seqname != transcript_id:
                if 'TRANSCRIPT_ID' in variant.attrs:
                    tx_id = variant.attrs['TRANSCRIPT_ID']
                    if tx_id == transcript_id:
                        seqname = tx_id
            if seqname not in variant_id_map:
                variant_id_map[seqname] = []
            variant_id_map[seqname].append(variant.id)
    if is_fusion:
        fusion_id = fusion_variant.id
        first_tx_id = fusion_variant.location.seqname
        first_gene_id = fusion_variant.attrs['GENE_ID']
        second_tx_id = fusion_variant.attrs['ACCEPTER_TRANSCRIPT_ID']
        second_gene_id = fusion_variant.attrs['ACCEPTER_GENE_ID']

        first_variant_ids = variant_id_map.get(first_gene_id, [])
        first_variant_ids += variant_id_map.get(first_tx_id, [])
        second_variant_ids = variant_id_map.get(second_gene_id, [])
        second_variant_ids += variant_id_map.get(second_tx_id, [])

        label = FusionVariantPeptideIdentifier(fusion_id, first_variant_ids,
            second_variant_ids, orf_id, index)
        return str(label)
    variant_ids = []
    for key, val in variant_id_map.items():
        if key == transcript_id or is_circ_rna:
            variant_ids += val
        else:
            variant_ids += [f"{key}-{x}" for x in val]
    if is_circ_rna:
        label = CircRNAVariantPeptideIdentifier(circ_rna_id, variant_ids, orf_id, index)
    else:
        label = BaseVariantPeptideIdentifier(transcript_id, variant_ids, orf_id, index)
    return str(label)

def parse_variant_peptide_id(label:str) -> List[VariantPeptideIdentifier]:
    """ Parse variant peptide identifier from peptide name.

    Examples:
        >>> peptide_label = 'ENST0001|SNV-100-A-T|1'
        >>> parse_variant_peptide_id(peptide_label)
        [<BaseVariantPeptideIdentifier>: 'ENST0001|SNV-100-A-T|1']

    """
    variant_ids = []
    for it in label.split(VARIANT_PEPTIDE_SOURCE_DELIMITER):
        x_id, *var_ids, index = it.split('|')

        orf_id = None
        if len(var_ids) > 0 and var_ids[0].startswith('ORF'):
            orf_id = var_ids.pop(0)

        if x_id.startswith('FUSION-'):
            first_variants:List[str] = []
            second_variants:List[str] = []
            for var_id in var_ids:
                which_gene, var_id = var_id.split('-', 1)
                if int(which_gene) == 1:
                    first_variants.append(var_id)
                elif int(which_gene) == 2:
                    second_variants.append(var_id)
                else:
                    raise ValueError('Variant is not valid')
            variant_id = FusionVariantPeptideIdentifier(x_id, first_variants,
                second_variants, orf_id, index)
        elif x_id.startswith('CIRC-') or x_id.startswith('CI-'):
            variant_id = CircRNAVariantPeptideIdentifier(x_id, var_ids, orf_id, index)
        else:
            variant_id = BaseVariantPeptideIdentifier(x_id, var_ids, orf_id, index)

        variant_ids.append(variant_id)
    return variant_ids


class VariantPeptideIdentifier(ABC):
    """ variant peptide identifer virtual class """
    @abstractmethod
    def __str__(self) -> str:
        """ str """

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>: '{str(self)}'"

class BaseVariantPeptideIdentifier(VariantPeptideIdentifier):
    """ Variant peptide identifier for output FASTA header """
    def __init__(self, transcript_id:str, variant_ids:List[str],
            orf_id:str=None, index:int=None):
        """ constructor """
        self.transcript_id = transcript_id
        self.variant_ids = variant_ids
        self.index = index
        self.orf_id = orf_id

    def __str__(self) -> str:
        """ str """
        x = [self.transcript_id]
        if self.orf_id:
            x.append(self.orf_id)
        x += self.variant_ids
        if self.index:
            x.append(str(self.index))
        return '|'.join(x)

    def is_alternative_splicing(self) -> bool:
        """ Whether this variant peptide has any alternative splicing events """
        alt_splice_types = ['SE', 'A5SS', 'A3SS', 'RI', 'MXE']
        return any(any(y in x for y in alt_splice_types) for x in self.variant_ids)

class CircRNAVariantPeptideIdentifier(VariantPeptideIdentifier):
    """ circRNA variant peptide identifier for output FASTA header """
    def __init__(self, circ_rna_id:str, variant_ids:List[str],
            orf_id:str=None, index:int=None):
        """ constructor """
        self.circ_rna_id = circ_rna_id
        self.variant_ids = variant_ids
        self.index = index
        self.orf_id = orf_id

    def __str__(self) -> str:
        """ str """
        x = [self.circ_rna_id]
        if self.orf_id:
            x.append(self.orf_id)
        x += self.variant_ids
        if self.index:
            x.append(str(self.index))
        return '|'.join(x)

class FusionVariantPeptideIdentifier(VariantPeptideIdentifier):
    """ Fusion variant peptide identifier """
    def __init__(self, fusion_id:str, first_variants:List[str],
            second_variants:List[str], orf_id:str=None, index:int=None):
        self.fusion_id = fusion_id
        self.first_variants = first_variants
        self.second_variants = second_variants
        self.index = index
        self.orf_id = orf_id

    def __str__(self) -> str:
        """ str """
        x = [self.fusion_id]
        if self.orf_id:
            x.append(self.orf_id)
        x += [f"1-{it}" for it in self.first_variants]
        x += [f"2-{it}" for it in self.second_variants]
        if self.index:
            x.append(str(self.index))
        return '|'.join(x)

    @property
    def first_tx_id(self):
        """ get first gene id """
        _,first,_ = self.fusion_id.split('-')
        return first.split(':')[0]

    @property
    def second_tx_id(self):
        """ get first gene id """
        _,_,second = self.fusion_id.split('-')
        return second.split(':')[0]
