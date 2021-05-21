"""Main entry point"""
import argparse
from moPepGen import vep


def main():
    """ Main entry point """
    parser = argparse.ArgumentParser(
        prog='moPopGen'
    )
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
        nargs='+',
        help='Path to VEP result txt file.',
        metavar=''
    )
    parser_vep.add_argument(
        '-g', '--genome-fasta',
        type=str,
        help='Path to the genome assembly FASTA file.',
        metavar=''
    )
    parser_vep.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.',
        metavar=''
    )
    parser_vep.add_argument(
        '-p', '--protein-fasta',
        type=str,
        help='Path to the protein FASTA file of the coding regions of the '
        'genome. Must come from ENSEMBL/GENCODE with the same version of the '
        'genome FASTA and GTF file.',
        metavar=''
    )
    parser_vep.add_argument(
        '-o', '--output-fasta',
        type=str,
        help='The output FASTA file name.',
        metavar=''
    )
    parser_vep.add_argument(
        '--cleavage-rule',
        type=str,
        help='The peptide cleavage rule for in silico digestion. Deafult to '
        'trypsin.',
        default='trypsin',
        metavar=''
    )
    parser_vep.add_argument(
        '--cleavage-exception',
        type=str,
        default=None,
        help='The peptide cleavage exception for in silico digestion. Default '
        'to None',
        metavar=''
    )
    parser_vep.add_argument(
        '--miscleavage',
        type=int,
        default=2,
        help='Number of miscleavages allowed for in silico digestion. Default '
        'to 2.',
        metavar=''
    )
    parser_vep.add_argument(
        '--min-molecular-weight',
        type=float,
        default=500.,
        help='The minimal molecular weight of peptides to report. Default to '
        '500.',
        metavar=''
    )
    parser_vep.set_defaults(func=vep.main)

    args = parser.parse_args()

    args.func(args)


if __name__ == '__main__':
    main()