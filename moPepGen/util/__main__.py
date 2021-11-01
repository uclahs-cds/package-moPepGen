"""Main entry point for moPepGen-util """
import sys
import argparse
from moPepGen import __version__, util


def main():
    """ Main entry point """

    parser = argparse.ArgumentParser(
        prog='moPopGen-util'
    )
    parser.add_argument(
        '-V', '--version',
        help='Show version number and exit.',
        action='store_true'
    )
    subparsers = parser.add_subparsers(
        dest='command'
    )

    util.add_subparser_brute_force(subparsers)
    util.add_subparser_downsample_reference(subparsers)
    util.add_subparser_validate_variant_callilng(subparsers)

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
