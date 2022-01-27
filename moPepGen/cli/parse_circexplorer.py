r""" `parseCIRCExplorer` takes the identified circRNA results from
[CIRCexplorer](https://circexplorer2.readthedocs.io/) and save as a
GVF file. The GVF file can be later used to call variant peptides using
[callVariant](call-variant.md). Noted that only known circRNA is supported (
\*_circular_known.txt) """
from __future__ import annotations
import argparse
from logging import warning
from typing import List, Dict
from pathlib import Path
from moPepGen import logger, circ, err
from moPepGen.parser import CIRCexplorerParser
from moPepGen.cli.common import add_args_input_path, add_args_reference, add_args_quiet, \
    add_args_source, add_args_output_path, print_start_message, \
    print_help_if_missing_args, load_references, generate_metadata, parse_range, validate_file_format


INPUT_FILE_FORMATS = ['.tsv', '.txt']
OUTPUT_FILE_FORMATS = ['.gvf']

# pylint: disable=W0212
def add_subparser_parse_circexplorer(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen parseCIRCexplorer """
    p = subparsers.add_parser(
        name='parseCIRCexplorer',
        help='Parse CIRCexplorer result',
        description='Parse CIRCexplorer result to a TSV format for moPepGen to'
        ' call variant peptides',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_args_input_path(
        parser=p, formats=INPUT_FILE_FORMATS,
        message="File path to CIRCexplorer's TSV output."
    )
    add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--circexplorer3',
        action='store_true',
        help='Using circRNA resutls from CIRCexplorer3'
    )
    p.add_argument(
        '--min-read-number',
        type=int,
        help='Minimal number of junction read counts.',
        default=1,
        metavar='<number>'
    )
    p.add_argument(
        '--min-fpb-circ',
        type=float,
        help='Minimal CRICscore value for CIRCexplorer3. Recommends to 1,'
        'defaults to None',
        default=None,
        metavar='<number>'
    )
    p.add_argument(
        '--min-circ-score',
        type=float,
        help='Minimal CIRCscore value for CIRCexplorer3. Recommends to 1,'
        'defaults to None',
        default=None,
        metavar='<number>'
    )
    p.add_argument(
        '--intron-start-range',
        type=str,
        help='The range of difference allowed between the intron start and'
        ' the reference position.',
        default='-2,0',
        metavar='<number>'
    )
    p.add_argument(
        '--intron-end-range',
        type=str,
        help='The range of difference allowed between the intron end and'
        ' the reference position.',
        default='-100,5',
        metavar='<number>'
    )
    add_args_source(p)
    add_args_reference(p, genome=False, proteome=False)
    add_args_quiet(p)
    p.set_defaults(func=parse_circexplorer)
    print_help_if_missing_args(p)
    return p

def parse_circexplorer(args:argparse.Namespace):
    """ Parse circexplorer known circRNA results. """
    input_path:Path = args.input_path
    output_path:Path = args.output_path
    validate_file_format(input_path, INPUT_FILE_FORMATS)
    validate_file_format(output_path, OUTPUT_FILE_FORMATS)

    intron_start_range = parse_range(args.intron_start_range)
    intron_end_range = parse_range(args.intron_end_range)

    print_start_message(args)

    _, anno, *_ = load_references(args, False, False)

    circ_records:Dict[str, List[circ.CircRNAModel]] = {}

    for record in CIRCexplorerParser.parse(input_path, args.circexplorer3):
        if not args.circexplorer3:
            if not record.is_valid(args.min_read_number):
                continue
        elif not record.is_valid(args.min_read_number, args.min_fbr_circ, \
                args.min_circ_score):
            continue
        try:
            circ_record = record.convert_to_circ_rna(anno, intron_start_range,
                intron_end_range)
        except err.ExonNotFoundError:
            err.warning(f"The CIRCexplorer record {record.name} from"
                f" transcript {record.isoform_name} contains an unknown exon."
                " Skipping it from parsing.")
            continue
        except err.IntronNotFoundError:
            err.warning(f"The CIRCexplorer record {record.name} from"
                f" transcript {record.isoform_name} contains an unknown intron."
                " Skipping it from parsing.")
            continue
        except:
            logger(f'Exception raised from record: {record.name}')
            raise
        gene_id = circ_record.gene_id
        if gene_id not in circ_records:
            circ_records[gene_id] = []
        circ_records[gene_id].append(circ_record)

    if not circ_records:
        if not args.quiet:
            warning('No variant record is saved.')
        return

    genes_rank = anno.get_genes_rank()
    ordered_keys = sorted(circ_records.keys(), key=lambda x:genes_rank[x])

    records = []
    for key in ordered_keys:
        val = circ_records[key]
        records.extend(val)

    metadata = generate_metadata(args)

    with open(output_path, 'w') as handle:
        circ.io.write(records, metadata, handle)

    if not args.quiet:
        logger("CircRNA records written to disk.")
