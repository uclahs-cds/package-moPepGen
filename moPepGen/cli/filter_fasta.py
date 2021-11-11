""" `filterFasta` takes the variant peptide sequence file (FASTA) and filters it
based on the gene expression data. A expresion table must be given as a CSV
or TSV. """
from __future__ import annotations
import argparse
from pathlib import Path
from typing import IO, Dict
from moPepGen.aa import VariantPeptidePool
from .common import add_args_reference, add_args_verbose, print_start_message,\
    print_help_if_missing_args, load_references, logger


# pylint: disable=W0212
def add_subparser_filter_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen splitDatabase """
    p = subparser.add_parser(
        name='filterFasta',
        help='Filter noncanonical peptides.',
        description='Filter noncanonical peptides according to gene expression'
        ' or gene biotypes.'
    )
    p.add_argument(
        '-i', '--input-fasta',
        type=Path,
        help='Input FASTA file, must be called by moPepGen.'
    )
    p.add_argument(
        '-o', '--output-fasta',
        type=Path,
        help='Output FASTA file.'
    )
    p.add_argument(
        '--exprs-table',
        type=Path,
        help="Path to the RNAseq quantification results."
    )
    p.add_argument(
        '--skip-lines',
        type=int,
        help='Number of lines to skip when reading the expression table.'
        'Defaults to 0',
        default=0
    )
    p.add_argument(
        '--delimiter',
        type=str,
        help='Delimiter of the expression table. Defaults to tab.',
        default='\t'
    )
    p.add_argument(
        '--tx-id-col',
        type=int,
        help="The index for transcript ID in the RNAseq quantification results."
        " Index is 1-based."
    )
    p.add_argument(
        '--quant-col',
        type=str,
        help='The column index number for quantification. Index is 1-based.'
    )
    p.add_argument(
        '--quant-cutoff',
        type=float,
        help='Quantification cutoff.'
    )
    p.add_argument(
        '--keep-all-coding',
        action='store_true',
        default=False,
        help='Keep all coding genes, regardless of their expression level.'
    )
    p.add_argument(
        '--keep-all-noncoding',
        action='store_true',
        default=False,
        help='Keep all noncoding genes, regardless of their expression level.'
    )

    add_args_reference(p, genome=False, proteome=False)
    add_args_verbose(p)
    print_help_if_missing_args(p)
    p.set_defaults(func=filter_fasta)
    return p

def filter_fasta(args:argparse.Namespace) -> None:
    """ Filter noncanonical peptide FASTA """
    print_start_message(args)

    _, anno, *_ = load_references(args, load_genome=False, \
        load_proteome=False, load_canonical_peptides=False)

    with open(args.input_fasta, 'rt') as handle:
        pool = VariantPeptidePool.load(handle)

    if args.verbose:
        logger('Peptide FASTA file loaded.')

    with open(args.exprs_table, 'rt') as handle:
        exprs = load_expression_table(
            handle=handle,
            tx_col=args.tx_id_col - 1,
            quant_col=args.quant_col - 1,
            skip=args.skip_lines,
            delim=args.delimiter
        )

    if args.verbose:
        logger('Gene expression table loaded.')

    filtered_pool = pool.filter(exprs, args.quant_cutoff, anno,
        args.keep_all_noncoding, args.keep_all_coding)

    filtered_pool.write(args.output_fasta)

    if args.verbose:
        logger('Filtered FASTA file saved.')

def load_expression_table(handle:IO, tx_col:int,quant_col:int, skip:int=0,
        delim:str='\t') -> Dict[str,float]:
    """ Load the gene expression quantification table """
    i = 0
    while i < skip:
        handle.readline()
        i += 1

    data = {}

    line:str
    for line in handle:
        fields = line.rstrip().split(delim)
        tx_id = fields[tx_col]
        quant = fields[quant_col]
        data[tx_id] = quant

    return data
