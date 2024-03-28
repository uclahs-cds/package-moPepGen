"""Main entry point"""
import sys
import argparse
from moPepGen import cli, __version__, constant, get_logger
from moPepGen.util.ResourcesMonitor import ResourcesMonitor
from moPepGen.cli.common import setup_loggers

CLI_MAIN_USAGE = f"{constant.PROG_NAME} [-h] [-V] <command> [options]"
CLI_MAIN_DESCRIPTION = """
-- Indexing
   generateIndex       Generate genome and proteome index files for moPepGen.
   updateIndex         Update moPepGen index.
   indexGVF            Generate index for GVF files.

-- Parsing
   parseVEP            Parse VEP output.
   parseREDItools      Parse REDItools annotated output.
   parseSTARFusion     Parse STAR-Fusion output.
   parseFusionCatcher  Parse FusionCatcher output.
   parseArriba         Parse Arriba output.
   parseRMATS          Parse rMATS output.
   parseCIRCexplorer   Parse CIRCexplorer known circRNA output.

-- Calling
   callVariant         Call non-canonical peptides from genomic variants.
   callNovelORF        Call non-canonical peptides from novel ORFs.
   callAltTranslation  Call non-canonital peptides with alternative translation
                       from coding transcripts.

-- Processing
   summarizeFasta      Summarize variant peptides.
   filterFasta         Filter noncanonical peptides.
   splitFasta          Split noncanonical peptides into separate databases.
   mergeFasta          Merge multiple variant peptide FASTA databases.
   encodeFasta         Encode variant peptide FASTA file header.
   decoyFasta          Generate decoy database FASTA file.
"""

def main():
    """ Main entry point """
    process_monitor = ResourcesMonitor()

    parser = argparse.ArgumentParser(
        prog=constant.PROG_NAME,
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
    cli.add_subparser_update_index(subparsers)
    cli.add_subparser_parse_vep(subparsers)
    cli.add_subparser_parse_reditools(subparsers)
    cli.add_subparser_parse_star_fusion(subparsers)
    cli.add_subparser_parse_fusion_catcher(subparsers)
    cli.add_subparser_parse_arriba(subparsers)
    cli.add_subparser_parse_rmats(subparsers)
    cli.add_subparser_parse_circexplorer(subparsers)
    cli.add_subparser_call_variant(subparsers)
    cli.add_subparser_call_novel_orf(subparsers)
    cli.add_subparser_call_alt_translation(subparsers)
    cli.add_subparser_split_fasta(subparsers)
    cli.add_subparser_filter_fasta(subparsers)
    cli.add_subparser_index_gvf(subparsers)
    cli.add_subparser_merge_fasta(subparsers)
    cli.add_subparser_encode_fasta(subparsers)
    cli.add_subparser_decoy_fasta(subparsers)
    cli.add_subparser_summarize_fasta(subparsers)

    # allowing values to start with -, such as -100,3
    # https://stackoverflow.com/a/21446783/11081630
    for i, arg in enumerate(sys.argv):
        if arg[0] == '-' and arg[1].isdigit():
            sys.argv[i] = ' ' + arg

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    debug_level = 0 if args.quiet else args.debug_level
    setup_loggers(debug_level)
    logger = get_logger()

    if args.version:
        print(f'moPepGen {__version__}', file=sys.stdout, flush=True)
        sys.exit()

    args.func(args)

    resources_usage = process_monitor.get_resource_usage()
    logger.info(resources_usage)

if __name__ == '__main__':
    main()
