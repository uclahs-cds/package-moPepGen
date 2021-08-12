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

    cli.add_subparser_generate_index(subparsers)
    cli.add_subparser_parse_vep(subparsers)
    cli.add_subparser_parse_reditools(subparsers)
    cli.add_subparser_parse_star_fusion(subparsers)
    cli.add_subparser_parse_fusion_catcher(subparsers)
    cli.add_subparser_parse_rmats(subparsers)
    cli.add_subparser_call_peptides(subparsers)
    cli.add_subparser_parse_circexplorer(subparsers)

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args.func(args)


if __name__ == '__main__':
    main()
