""" VEP2VariantPeptides module """
import argparse
from Bio import SeqIO
from moPepGen.io import VepIO, GtfIO
from .VEP2VariantPeptides import VEP2VariantPeptides


def main(args) -> None:
    """ Main entry point for the VEP2VariantPeptides module. """
    peptides = VEP2
    peptides = VEP2VariantPeptides.setup(
        vep_path=args.vep_path,
        gtf_path=args.gtf_path,
        genome_path=args.genome_path,
        proteome_path=args.proteome_path
    ) \
        .call_canonical_peptides() \
        .call_variant_peptides()
    
    with open(args.output_fasta, 'r') as handle:
        for seq in peptides:
            SeqIO.write(seq, handle, 'fasta')
