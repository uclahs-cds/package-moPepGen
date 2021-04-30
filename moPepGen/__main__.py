"""Main entry point"""
import argparse
from moPepGen.vep import main as vep


def main():
    """ Main entry point """
    parser = argparse.ArgumentParser(prog='moPopGen')
    subparsers = parser.add_subparsers()

    ## VEP
    parser_vep = subparsers.add_parser(
        name='vep',
        help='Generate variant peptides from VEP output.',
        description="Generate variant peptides from VEP output. The same"
        "ENSEMBL/GENCODE version must be used for the genome assembly FASTA, "
        "protein sequence FASTA, and GTF files."
    )

    parser_vep.add_argument(
        '-v', '--vep-txt',
        type=str,
        help='Path to VEP result txt file.'
    )
    parser_vep.add_argument(
        '-g', '--genome-fasta',
        type=str,
        help='Path to the genome assembly FASTA file.'
    )
    parser_vep.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.'
    )
    parser_vep.add_argument(
        '-p', '--protein-fasta',
        type=str,
        help='Path to the protein FASTA file of the coding regions of the '
        'genome. Must come from ENSEMBL/GENCODE with the same version of the '
        'genome FASTA and GTF file.'
    )
    parser_vep.set_defaults(func=vep)

    args = parser.parse_args()

    args.func(args)


if __name__ == '__main__':
    main()