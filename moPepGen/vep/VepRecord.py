""" This module defines the class for a single VEP record
"""
from typing import List, Tuple
import re
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, dna


class VEPRecord():
    """ A VEPRecord object holds the an entry from the VEP output. The VEP
    output is defined at https://uswest.ensembl.org/info/docs/tools/vep/
    vep_formats.html#output
    
    Attributes:
        uploaded_variation (str): as chromosome_start_alleles
        location (str): in standard coordinate format (chr:start or
            chr:start-end)
        allele (str): the variant allele used to calculate the consequence
        gene (str): Ensembl stable ID of affected gene
        feature (str): Ensembl stable ID of feature
        feature_type (str): type of feature. Currently one of Transcript,
            RegulatoryFeature, MotifFeature.
        consequence (List[str]): consequence type of this variant.
            See: https://uswest.ensembl.org/info/genome/variation/prediction/
            predicted_data.html#consequences
        cdna_position (str): relative position of base pair in cDNA sequence
        cds_position (str): relative position of base pair in coding sequence
        protein_position (str): relative position of amino acid in protein
        amino_acids (Tuple[str]): only given if the variant affects the
            protein-coding sequence
        codons (Tuple[str]) the alternative codons with the variant base in
            upper case
        existing_variation (str) known identifier of existing variant
        extra (dict): this column contains extra information.
    """
    def __init__(
            self, uploaded_variation: str, location: str, allele: str,
            gene: str, feature: str, feature_type:str,
            consequences: List[str], cdna_position: int, cds_position: int,
            protein_position: int, amino_acids: Tuple[str, str],
            codons: Tuple[str, str], existing_variation: str, extra: dict):
        """ Construct a VEPRecord object. """
        self.uploaded_variation = uploaded_variation
        self.location = location
        self.allele = allele
        self.gene = gene
        self.feature = feature
        self.feature_type=feature_type
        self.consequences = consequences
        self.cdna_position = cdna_position
        self.cds_position = cds_position
        self.protein_position = protein_position
        self.amino_acids = amino_acids
        self.codons = codons
        self.existing_variation = existing_variation
        self.extra = extra
    
    def __repr__(self)->str:
        """Return representation of the VEP record."""
        consequences = '|'.join(self.consequences)
        return f"< {self.feature}, {consequences}, {self.location} >"
    
    def convert_to_variant_record(self, seq:dna.DNASeqRecord
            ) -> seqvar.VariantRecord:
        """ """
        alt_position = self.cdna_position.split('-')
        alt_start = int(alt_position[0]) - 1
        
        codon_ref, codon_alt = self.codons

        if codon_ref == '-':
            alt_end = alt_start + 1
            ref = seq.seq[alt_start:alt_end]
            alt = ref + codon_alt
            alt_end = alt_start + len(ref)
        elif codon_alt == '-':
            alt_start -= 1
            alt_end = alt_start + len(codon_ref) + 1
            ref = seq.seq[alt_start:alt_end]
            alt = seq.seq[alt_start:alt_start+1]
        else:
            pattern = re.compile('[ATCG]+')
            match_ref = pattern.search(codon_ref)
            match_alt = pattern.search(codon_alt)
            if match_ref is not None:
                if match_alt is not None:
                    ref = match_ref.group()
                    alt_end = alt_start + len(ref)
                    alt = match_alt.group()
                else:
                    ref = match_ref.group()
                    alt_start -= 1
                    ref = seq.seq[alt_start] + ref
                    alt_end = alt_start + len(ref)
                    alt = seq.seq[alt_start:alt_end]
            elif match_alt is not None:
                alt = match_alt.group()
                # alt_start -= 1
                alt_end = alt_start + 1
                ref = seq.seq[alt_start:alt_end]
                alt = ref + alt
            else:
                raise ValueError('No alteration found in this VEP record')

        type = 'SNV' if len(ref) == 1 and len(alt) == 1 else 'INDEL'
        _id = f'{type}-{alt_start}-{ref}-{alt}'

        try:
            return seqvar.VariantRecord(
                location=FeatureLocation(
                    seqname=self.feature,
                    start=alt_start,
                    end=alt_end
                ),
                ref=ref,
                alt=alt,
                type=type,
                id=_id
            )
        except ValueError as e:
            raise ValueError(e.args[0] + f' [{self.feature}]')
        