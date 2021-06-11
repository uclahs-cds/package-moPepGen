"""Main entry point"""
import sys
import argparse
from moPepGen import cli


def main():
    """ Main entry point """
    parser = argparse.ArgumentParser(
        prog='moPopGen'
    )
    subparsers = parser.add_subparsers()

    ## generate index
    parser_index = subparsers.add_parser(
        name='generateIndex',
        help='Generate genome and proteome index files for moPepGen',
        description='Generate genome and proteome index files for moPepGen'
        'parsers and peptide caller.'
    )
    parser_index.add_argument(
        '-g', '--genome-fasta',
        type=str,
        help='Path to the genome assembly FASTA file.',
        metavar='',
        required=True
    )
    parser_index.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.',
        metavar='',
        required=True
    )
    parser_index.add_argument(
        '-p', '--proteome-fasta',
        type=str,
        help='Path to the translated protein sequence FASTA file. Must come'
        'from ENSEMBL/GENCODE with the same version of the genome FASTA.',
        metavar='',
        required=True
    )
    parser_index.add_argument(
        '-c', '--cleavage-rule',
        type=str,
        help='Cleavage rule. Defaults to trypsin.',
        default='trypsin',
        metavar=''
    )
    parser_index.add_argument(
        '-m', '--miscleavage',
        type=int,
        help='Number of cleavages to allow. Defaults to 2.',
        metavar='',
        default=2
    )
    parser_index.add_argument(
        '-w', '--min-mw',
        type=float,
        help='The minimal molecular weight of the non-canonical peptides.'
        'Defaults to 500',
        default=500.,
        metavar=''
    )
    parser_index.add_argument(
        '-o', '--output-dir',
        type=str,
        help='Ouput directory for index files.',
        metavar='',
        dest='output_dir',
        required=True
    )
    parser_index.add_argument(
        '--verbose',
        type=bool,
        help='Print out logging messages. Defaults to True',
        default=True
    )
    parser_index.set_defaults(func=cli.generate_index)

    ## parseVEP
    parser_parse_vep = subparsers.add_parser(
        name='parseVEP',
        help='Parse VEP output for moPepGen to call variant peptides.',
        description="Parse VEP output tsv to a BED-like format of variant"
        "information for moPepGen to call variant peptides. The genome"
        "assembly FASTA and annotation GTF must come from the same"
        "GENCODE/ENSEMBL version, and must the consistent with the VEP output."
    )

    parser_parse_vep.add_argument(
        '-v', '--vep-txt',
        type=str,
        nargs='+',
        help='Path to VEP result txt file.',
        metavar='',
        required=True
    )
    parser_parse_vep.add_argument(
        '-g', '--genome-fasta',
        type=str,
        help='Path to the genome assembly FASTA file.',
        metavar=''
    )
    parser_parse_vep.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.',
        metavar=''
    )
    parser_parse_vep.add_argument(
        '--index-dir',
        type=str,
        help='Path to the directory of index files generated by moPepGen'
        'generateIndex. If given, --geome-fasta and --anntotation-gtf will be'
        'ignored.',
        metavar='',
        nargs='?',
        default=None
    )
    parser_parse_vep.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar='',
        required=True
    )
    parser_parse_vep.add_argument(
        '--verbose',
        type=bool,
        help='Print out logging messages. Defaults to True',
        default=True
    )

    parser_parse_vep.set_defaults(func=cli.parse_vep)

    ## parseREDItools
    parser_parse_reditools = subparsers.add_parser(
        name='parseREDItools',
        help='Parse REDItools result for moPepGen to call variant peptides.',
        description='Parse the REDItools result to a BED-like format of'
        'variant information for moPepGen to call variant peptides. The genome'
    )
    parser_parse_reditools.add_argument(
        '-t', '--reditools-table',
        type=str,
        help='Path to the REDItools output table.',
        metavar='',
        required=True
    )
    parser_parse_reditools.add_argument(
        '--transcript-id-column',
        type=int,
        help='The column index for transcript ID. If your REDItools table does'
        'not contains it, use the AnnotateTable.py from the REDItools'
        'package. Defaults to 16',
        default=16,
        metavar=''
    )
    parser_parse_reditools.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Ignored if --index-dir is'
        'specified',
        metavar=''
    )
    parser_parse_reditools.add_argument(
        '--index-dir',
        type=str,
        help='Path to the directory of index files generated by moPepGen'
        'generateIndex. If given, --geome-fasta and --anntotation-gtf will be'
        'ignored.',
        metavar='',
        nargs='?',
        default=None
    )
    parser_parse_reditools.add_argument(
        '-o', '--output-prefix',
        type=str,
        help='Prefix to the output filename.',
        metavar='',
        required=True
    )
    parser_parse_reditools.add_argument(
        '--verbose',
        type=bool,
        help='Print out logging messages. Defaults to True',
        default=True
    )

    parser_parse_reditools.set_defaults(func=cli.parse_reditools)


    ## callPeptides
    parser_call_peptide = subparsers.add_parser(
        name='callPeptide',
        help='Call non-canonical peptides from genomic variants.',
        description="Genomic variant data must be generated by one of the"
        "moPepGen parser. See moPepGen --help"
    )

    parser_call_peptide.add_argument(
        '-i', '--input-variant',
        type=str,
        nargs='+',
        help='Path to input variant files. Must be generated by any of the'
        'moPepGen parser. This can be multiple.',
        metavar='',
        required=True
    )
    parser_call_peptide.add_argument(
        '-g', '--genome-fasta',
        type=str,
        help='Path to the genome assembly FASTA file.',
        metavar=''
    )
    parser_call_peptide.add_argument(
        '-a', '--annotation-gtf',
        type=str,
        help='Path to the annotation GTF file. Must come from ENSEMBL/GENCODE'
        ' with the same version of the genome and protein FASTA.',
        metavar=''
    )
    parser_call_peptide.add_argument(
        '-p', '--proteome-fasta',
        type=str,
        help='Path to the translated protein sequence FASTA file. Must come'
        'from ENSEMBL/GENCODE with the same version of the genome FASTA.',
        metavar=''
    )
    parser_call_peptide.add_argument(
        '--index-dir',
        type=str,
        help='Path to the directory of index files generated by moPepGen'
        'generateIndex. If given, --geome-fasta, --proteome-fasta and'
        '--anntotation-gtf will be ignored.',
        metavar='',
        nargs='?',
        default=None
    )
    parser_call_peptide.add_argument(
        '-o', '--output-fasta',
        type=str,
        help='Filename for the output FASTA.',
        metavar='',
        required=True
    )
    parser_call_peptide.add_argument(
        '-c', '--cleavage-rule',
        type=str,
        help='Cleavage rule. Defaults to trypsin. Defaults to trypsin',
        default='trypsin',
        metavar=''
    )
    parser_call_peptide.add_argument(
        '-m', '--miscleavage',
        type=int,
        help='Number of cleavages to allow. Defaults to 2',
        default=2,
        metavar=''
    )
    parser_call_peptide.add_argument(
        '-w', '--min-mw',
        type=float,
        help='The minimal molecular weight of the non-canonical peptides.'
        'Defaults to 500',
        default=500.,
        metavar=''
    )
    parser_call_peptide.add_argument(
        '-v', '--verbose',
        type=bool,
        help='Verbose',
        metavar='',
        default=True
    )

    parser_call_peptide.set_defaults(func=cli.call_variant_peptide)

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args.func(args)


if __name__ == '__main__':
    main()
