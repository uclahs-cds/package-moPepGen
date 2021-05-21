""" VEP2VariantPeptides module """
import argparse
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord
from moPepGen.vep.VEP2VariantPeptides import VEP2VariantPeptides


def main(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP2VariantPeptides module. """
    adapter = VEP2VariantPeptides.set_up(
        vep_path=args.vep_txt,
        gtf_path=args.annotation_gtf,
        genome_path=args.genome_fasta,
        proteome_path=args.protein_fasta
    )
    adapter.call_variant_peptides(
        rule=args.cleavage_rule,
        exception=args.cleavage_exception,
        miscleavage=args.miscleavage,
        min_mw=args.min_molecular_weight
    )
    adapter.write_peptides(args.output_fasta)
