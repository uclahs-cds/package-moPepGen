"""Main entry point"""
import sys
import argparse
from moPepGen import cli, __version__

CLI_MAIN_USAGE = "moPopGen [-h] [-V] <command> [options]"
CLI_MAIN_DESCRIPTION = """
-- Indexing
   generateIndex       Generate genome and proteome index files for moPepGen.
   indexGVF            Generate index for GVF files.

-- Parsing
   parseVEP            Parse VEP output.
   parseREDItools      Parse REDItools annotated output.
   parseSTARFusion     Parse STAR-Fusion output.
   parseFusionCatcher  Parse FusionCatcher output.
   parseArriba        Parse Arriba output.
   parseRMATS          Parse rMATS output.
   parseCIRCexplorer   Parse CIRCexplorer known circRNA output.

-- Calling
   callVariant         Call non-canonical peptides from genomic variants.
   callNoncoding       Call non-canonical peptides from noncoding transcripts.

-- Processing
   filterFasta         Filter noncanonical peptides.
   splitDatabase       Split noncanonical peptides into separate databases.
"""

def main():
    """ Main entry point """

    parser = argparse.ArgumentParser(
        prog='moPopGen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=CLI_MAIN_USAGE
    )
    parser.add_argument(
        '-V', '--version',
        help='Show version number and exit.',
        action='store_true'
    )
    subparsers = parser.add_subparsers(
        dest='command',
        help=argparse.SUPPRESS,
        description=CLI_MAIN_DESCRIPTION
    )

    cli.add_subparser_generate_index(subparsers)
    cli.add_subparser_parse_vep(subparsers)
    cli.add_subparser_parse_reditools(subparsers)
    cli.add_subparser_parse_star_fusion(subparsers)
    cli.add_subparser_parse_fusion_catcher(subparsers)
    cli.add_subparser_parse_arriba(subparsers)
    cli.add_subparser_parse_rmats(subparsers)
    cli.add_subparser_parse_circexplorer(subparsers)
    cli.add_subparser_call_variant(subparsers)
    cli.add_subparser_call_noncoding(subparsers)
    cli.add_subparser_split_database(subparsers)
    cli.add_subparser_filter_fasta(subparsers)
    cli.add_subparser_index_gvf(subparsers)

    # allowing values to start with -, such as -100,3
    # https://stackoverflow.com/a/21446783/11081630
    for i, arg in enumerate(sys.argv):
        if arg[0] == '-' and arg[1].isdigit():
            sys.argv[i] = ' ' + arg

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.version:
        print(f'moPepGen {__version__}', file=sys.stdout, flush=True)
        sys.exit()

    args.func(args)


if __name__ == '__main__':
    main()
