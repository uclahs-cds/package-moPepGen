""" Module for VariantPeptideIdentifier """
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Dict, List, TYPE_CHECKING, Set
from moPepGen import VARIANT_PEPTIDE_SOURCE_DELIMITER
from moPepGen.constant import VariantPrefix

if TYPE_CHECKING:
    from moPepGen.seqvar import VariantRecord

def create_variant_peptide_id(transcript_id:str, variants:List[VariantRecord],
        orf_id:str=None, index:int=None, gene_id:str=None) -> str:
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
            if variant.is_merged_mnv():
                variant_id_map[seqname].extend(variant.attrs['INDIVIDUAL_VARIANT_IDS'])
            else:
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
            second_variant_ids, [], orf_id, index)
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
        label = BaseVariantPeptideIdentifier(transcript_id, variant_ids, orf_id, index, gene_id)
    return str(label)

def _set_identifier_type(t, v):
    """ Set identifier type """
    if t is not None:
        raise ValueError(f"VariantType is not None: {t}")
    return v

def parse_variant_peptide_id(label:str, coding_txs:Set[str]) -> List[VariantPeptideIdentifier]:
    """ Parse variant peptide identifier from peptide name.

    Examples:
        >>> peptide_label = 'ENST0001|SNV-100-A-T|1'
        >>> parse_variant_peptide_id(peptide_label)
        [<BaseVariantPeptideIdentifier>: 'ENST0001|SNV-100-A-T|1']

    """
    # FASTA header examples:
    # ENST0001|SNV-50-A-T|1
    # ENST0001|SNV-50-A-T|INDEL-55-CC-C|2
    # ENST0001|SNV-50-A-T|ORF2|1
    # FASTA headers have 5 fields below separate by | :
    # - (Required) Transcript backbone ID (tx ID, fusion ID, or circRNA ID)
    # - (Optioanl) Gene ID, required for novel ORF peptides without additional variants.
    # - (Optional) Variant IDs separated by |
    # - (Optional) ORF ID, required for novel ORF peptides.
    # - (Required) peptide index.
    variant_ids = []
    for it in label.split(VARIANT_PEPTIDE_SOURCE_DELIMITER):
        fields = it.split('|')
        gene_id = None
        backbone_id = None
        var_ids:Dict[int, List[str]] = {}
        alt_ids = []
        orf_id = None

        try:
            index = int(fields[-1])
            fields.pop()
        except ValueError:
            index = None

        IdentifierType:VariantPeptideIdentifier = None

        # Step 1. Iterate over the fields in a FASTA header separated by | and
        # assign them into the 5 variables.
        for i, field in enumerate(fields):
            # The first field is always the backbone ID, and if its prefix is
            # FUSION-, CIRC-, or CI-, this must be a Fusion or circRNA entry.
            if field.startswith(str(VariantPrefix.FUSION)):
                if i != 0:
                    raise ValueError(f"FASTA header isn't valid: {it}")
                backbone_id = field
                IdentifierType = _set_identifier_type(
                    IdentifierType,
                    FusionVariantPeptideIdentifier
                )

            elif field.startswith(str(VariantPrefix.CI)) \
                    or field.startswith(str(VariantPrefix.CIRC)):
                if i != 0:
                    raise ValueError(f"FASTA header isn't valid: {it}")
                backbone_id = field
                IdentifierType = _set_identifier_type(
                    IdentifierType,
                    CircRNAVariantPeptideIdentifier
                )

            elif field.startswith('ORF'):
                orf_id = field

            elif IdentifierType is FusionVariantPeptideIdentifier:
                if field.startswith('1-') or field.startswith('2-'):
                    which_gene, var_id = field.split('-', 1)
                    which_gene = int(which_gene)
                else:
                    # SECT or W2F
                    var_id = field
                    which_gene = 0
                if which_gene in var_ids:
                    var_ids[which_gene].append(var_id)
                else:
                    var_ids[which_gene] = [var_id]

            elif any(field.startswith(str(t)) for t in VariantPrefix.alt_translation()):
                alt_ids.append(field)

            elif any(field.startswith(str(t)) for t in VariantPrefix.ctbv()):
                if 1 in var_ids:
                    var_ids[1].append(field)
                else:
                    var_ids[1] = [field]

        # Step 2. Infer the FASTA header type. If it's not already interpreted in
        # step 1 (either Fusion or circRNA), it must be either a base variant or
        # novel ORF peptide.
        if IdentifierType is None:
            # `NovelORFPeptideIdentifier` is for novel ORF peptides without any
            # additional variants. While `BaseVariantPeptideIdentifer` always
            # have at least one variant.
            if not var_ids and orf_id is not None:
                backbone_id = fields[0]
                if len(fields) > 1:
                    gene_id = fields[1]
                IdentifierType = _set_identifier_type(
                    IdentifierType,
                    NovelORFPeptideIdentifier
                )
            else:
                backbone_id = fields[0]
                IdentifierType = _set_identifier_type(
                    IdentifierType,
                    BaseVariantPeptideIdentifier
                )

        # Step 3. Construct the object based on the FASTA header type.
        if IdentifierType is NovelORFPeptideIdentifier:
            variant_id = NovelORFPeptideIdentifier(
                transcript_id=backbone_id,
                gene_id=gene_id,
                codon_reassigns=alt_ids,
                orf_id=orf_id,
                index=index,
                is_protein_coding=(backbone_id in coding_txs)
            )
        elif IdentifierType is FusionVariantPeptideIdentifier:
            variant_id = FusionVariantPeptideIdentifier(
                fusion_id=backbone_id,
                first_variants=var_ids.get(1, []),
                second_variants=var_ids.get(2, []),
                peptide_variants=var_ids.get(0, []),
                orf_id=orf_id,
                index=index
            )
        elif IdentifierType is CircRNAVariantPeptideIdentifier:
            variant_id = CircRNAVariantPeptideIdentifier(
                circ_rna_id=backbone_id,
                variant_ids=sum(var_ids.values(), []) + alt_ids,
                orf_id=orf_id,
                index=index
            )
        elif IdentifierType is BaseVariantPeptideIdentifier:
            variant_id = BaseVariantPeptideIdentifier(
                transcript_id=backbone_id,
                variant_ids=sum(var_ids.values(), []) + alt_ids,
                orf_id=orf_id,
                index=index,
                gene_id=gene_id
            )

        variant_ids.append(variant_id)
    return variant_ids


class VariantPeptideIdentifier(ABC):
    """ variant peptide identifier virtual class """
    @abstractmethod
    def __str__(self) -> str:
        """ str """

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>: '{str(self)}'"

class BaseVariantPeptideIdentifier(VariantPeptideIdentifier):
    """ Variant peptide identifier for output FASTA header """
    def __init__(self, transcript_id:str, variant_ids:List[str],
            orf_id:str=None, index:int=None, gene_id:str=None):
        """ constructor """
        self.transcript_id = transcript_id
        self.variant_ids = variant_ids
        self.index = index
        self.orf_id = orf_id
        self.gene_id = gene_id

    def __str__(self) -> str:
        """ str """
        x = [self.transcript_id]
        if self.gene_id:
            x.append(self.gene_id)
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
            second_variants:List[str], peptide_variants:List[str],
            orf_id:str=None, index:int=None):
        self.fusion_id = fusion_id
        self.first_variants = first_variants
        self.second_variants = second_variants
        self.peptide_variants = peptide_variants
        self.index = index
        self.orf_id = orf_id

    def __str__(self) -> str:
        """ str """
        x = [self.fusion_id]
        if self.orf_id:
            x.append(self.orf_id)
        x += [f"1-{it}" for it in self.first_variants]
        x += [f"2-{it}" for it in self.second_variants]
        x += self.peptide_variants
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

class NovelORFPeptideIdentifier(VariantPeptideIdentifier):
    """ Novel ORF peptide identifier """
    def __init__(self, transcript_id:str, gene_id:str, codon_reassigns:List[str],
            is_protein_coding:bool, orf_id:str=None, index:int=None):
        """ constructor """
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.codon_reassigns = codon_reassigns
        self.orf_id = orf_id
        self.index = index
        self.is_protein_coding = is_protein_coding

    def __str__(self) -> str:
        """ str """
        fields = [self.transcript_id]
        if self.gene_id:
            fields.append(self.gene_id)
        fields += self.codon_reassigns
        if self.orf_id:
            fields.append(self.orf_id)
        if self.index:
            fields.append(str(self.index))
        return '|'.join(fields)
