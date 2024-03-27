"""Main entry point for moPepGen-util """
import sys
import argparse
from moPepGen import __version__, util
from moPepGen.cli.common import setup_loggers


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

    util.brute_force.parse_args(subparsers)
    util.brute_force_novel_orf.parse_args(subparsers)
    util.downsample_reference.parse_args(subparsers)
    util.validate_variant_calling.parse_args(subparsers)
    util.fuzz_test.parse_args(subparsers)
    util.extract_gvf.parse_args(subparsers)
    util.validate_novel_orf_calling.parse_args(subparsers)
    util.add_fuzz_test_log.parse_args(subparsers)

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    setup_loggers(args.debug_level)

    if args.version:
        print(f'moPepGen {__version__}', file=sys.stdout, flush=True)
        sys.exit()

    args.func(args)


if __name__ == '__main__':
    main()
