""" Module for common functions for the util module """
from pathlib import Path
from typing import Tuple
from moPepGen.gtf import GenomicAnnotation
from moPepGen.dna import DNASeqDict
from moPepGen.aa import AminoAcidSeqDict


def load_references(path_anno:Path, path_genome:Path, path_proteome:Path
        ) -> Tuple[GenomicAnnotation, DNASeqDict, AminoAcidSeqDict]:
    """ Load reference files """
    anno = GenomicAnnotation()
    anno.dump_gtf(path_anno)

    genome = DNASeqDict()
    genome.dump_fasta(path_genome)

    proteome = AminoAcidSeqDict()
    proteome.dump_fasta(path_proteome)

    anno.check_protein_coding(proteome, invalid_protein_as_noncoding=True)

    return anno, genome, proteome
