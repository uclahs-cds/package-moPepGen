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
    variant_ids:Dict[str,List[VariantRecord]] = {}
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
            gene_id = variant.location.seqname
            if gene_id not in variant_ids:
                variant_ids[gene_id] = []
            variant_ids[gene_id].append(variant.id)
    if is_fusion:
        fusion_id = fusion_variant.id
        first_gene_id = fusion_variant.location.seqname
        second_gene_id = fusion_variant.attrs['ACCEPTER_GENE_ID']
        label = FusionVariantPeptideIdentifer(fusion_id,
            variant_ids[first_gene_id], variant_ids[second_gene_id],
            index)
        return str(label)
    gene_ids = list(variant_ids.keys())
    if len(gene_ids) > 1:
        raise ValueError('Variants should all have the same gene ID')
    if len(gene_ids) == 0:
        variant_ids = []
    else:
        variant_ids = variant_ids[gene_ids[0]]
    if is_circ_rna:
        label = CircRNAVariantPeptideIdentifier(circ_rna_id, variant_ids, orf_id, index)
    else:
        label = BaseVariantPeptideIdentifier(transcript_id, variant_ids, orf_id, index)
    return str(label)

def parse_variant_peptide_id(label:str) -> List[VariantPeptideIdentifier]:
    """ Parse variant peptide info from label """
    variant_ids = []
    for it in label.split(VARIANT_PEPTIDE_SOURCE_DELIMITER):
        x_id, *var_ids, index = it.split('|')

        orf_id = None
        if len(var_ids) > 0 and var_ids[0].startswith('ORF'):
            orf_id = var_ids.pop(0)

        if x_id.startswith('FUSION'):
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
            variant_id = FusionVariantPeptideIdentifer(x_id, first_variants,
                second_variants, orf_id, index)
        elif '-circRNA-' in x_id:
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

class FusionVariantPeptideIdentifer(VariantPeptideIdentifier):
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
    def first_gene_id(self):
        """ get first gene id """
        _,first,_ = self.fusion_id.split('-')
        return first.split(':')[0]

    @property
    def second_gene_id(self):
        """ get first gene id """
        _,_,second = self.fusion_id.split('-')
        return second.split(':')[0]
