""" The module defines the main logic that calls variant peptides from the VEP
output.
"""
from __future__ import annotations
from typing import Dict, Set
from moPepGen.io import VepIO
from moPepGen.DNASeqDict import DNASeqDict
from moPepGen.AminoAcidSeqDict import AminoAcidSeqDict, AminoAcidSeqRecord
from moPepGen.TranscriptVEPDict import TranscriptVEPDict
from moPepGen.TranscriptGTFDict import TranscriptGTFDict


class VEP2VariantPeptides():
    """ Defines the variant peptide caller from VEP output.

    Attributes:
        vep (TranscriptVEPDict): A dict-like object that contains the variation
            effects predicted by VEP. Keys are transcript IDs and values are
            list of VEPRecord.
        annotation (TranscriptGTFDict): A dict-like object that contains the
            annotation parsed from a GTF file. Keys are transcript IDs and
            values are liks of GTFRecord.
        genome (Dict[str, [SeqRecord]])
        proteome (Dict[str, [SeqRecord]])
        canonical_peptides (set(SeqRecord))
    """
    def __init__(self,
            vep: TranscriptVEPDict,
            annotation: TranscriptGTFDict,
            genome: DNASeqDict,
            proteome: AminoAcidSeqDict,
            canonical_peptides: Set[Seq]=None):
        """ Create a VEP2VariantPeptides object.

        Args:
            vep (TranscriptVEPDict): A dict-like object that contains the
                variation effects predicted by VEP. Keys are transcript IDs
                and values are list of VEPRecord.
            annotation (TranscriptGTFDict): A dict-like object that contains
                the annotation parsed from a GTF file. Keys are transcript IDs
                and values are liks of GTFRecord.
            genome (Dict[str, [SeqRecord]])
            proteome (Dict[str, [SeqRecord]])
            canonical_peptides (set(SeqRecord))
        """
        self.vep = vep
        self.annotation = annotation
        self.genome = genome
        self.proteome = proteome
        self.canonical_peptides = canonical_peptides
    
    def dump_data_files(self, vep_path:str, gtf_path:str, genome_path:str,
            proteome_path:str)->None:
        """ Load VEP, GTF, and FASTA files from disk.
        """
        self.vep.dump_vep()
        self.annotation.dump_gtf()
        self.genome.dump_fasta()
        self.proteome.dump_fasta()

    @staticmethod
    def set_up(vep_path:str, gtf_path:str, genome_path:str,
            proteome_path:str)->VEP2VariantPeptides:
        """ This should be the main starting point of running the
        VEP2VariantPeptides module.

        Args:
            vep_path (str): Path to the VEP file.
            gtf_path (str): Path to the annotation GTF.
            genome_path (str): Path to the genome assembly FASTA.
            proteome_path (str): path to the translated FASTA of the coding
                transcripts of the genome assembly.
        
        Example:
            adapter = VEP2VariantPaptides.set_up(
                vep_path='path/to/vep.txt',
                gtf_path='path/to/annotation.gtf',
                genome_path='path/to/genome_assembly.fasta',
                proteome_path='path/to/proteins.fasta'
            )
        """
        adapter = VEP2VariantPeptides(
            vep=None,
            annotation=None,
            genome=None,
            proteome=None
        )
        adapter.dump_data_files(
            vep_path=vep_path,
            gtf_path=gtf_path,
            genome_path=genome_path,
            proteome_path=proteome_path
        )
        return adapter


    def call_canonical_peptides(self, enzyme:str='trypsin'):
        """ Performs in silico digestion on each protein sequence from the
        proteome and create a set of all unique peptides.

        Args:
            enzyme (str): The enzyme for in silico digestion. Defaults to
                trypsin.
        
        Returns:
            The original instance of VEP2VariantPeptides with the canonical
            peptides stored in the canonical_peptides attribute.
        """
        return self

    def call_variant_peptides(self, enzeym:str='trypsin'
            )->set[AminoAcidSeqRecord]:
        """ Performs in silico digestion on each variated peptide and create
        a set of unique peptides that do not overlap with the canonical
        peptide pool.

        Args:
            enzyme (str): The enzyme for in silico digestion. Defaults to
                trypsin.
        """
        pass

