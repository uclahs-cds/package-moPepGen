""" `filterFasta` takes the variant peptide sequence file (FASTA) and filters it
based on the gene expression data. A expresion table must be given as a CSV
or TSV. """
from __future__ import annotations
import argparse
from pathlib import Path
import pickle
from typing import IO, Dict, List
from moPepGen.aa import VariantPeptidePool
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation
from .common import add_args_reference, add_args_quiet, print_start_message,\
    print_help_if_missing_args, logger


# pylint: disable=W0212
def add_subparser_filter_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen splitDatabase """
    p = subparser.add_parser(
        name='filterFasta',
        help='Filter noncanonical peptides.',
        description='Filter noncanonical peptides according to gene expression'
        ' or gene biotypes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-i', '--input-fasta',
        type=Path,
        help='Input FASTA file, must be called by moPepGen.',
        metavar='<file>'
    )
    p.add_argument(
        '-o', '--output-fasta',
        type=Path,
        help='Output FASTA file.',
        metavar='<file>'
    )
    p.add_argument(
        '--exprs-table',
        type=Path,
        help="Path to the RNAseq quantification results.",
        metavar='<file>'
    )
    p.add_argument(
        '--skip-lines',
        type=int,
        help='Number of lines to skip when reading the expression table.'
        'Defaults to 0',
        metavar='<value>',
        default=0
    )
    p.add_argument(
        '--delimiter',
        type=str,
        help='Delimiter of the expression table. Defaults to tab.',
        metavar='<value>',
        default='\t'
    )
    p.add_argument(
        '--tx-id-col',
        type=str,
        help="The index for transcript ID in the RNAseq quantification results."
        " Index is 1-based.",
        metavar='<number>'
    )
    p.add_argument(
        '--quant-col',
        type=str,
        help='The column index number for quantification. Index is 1-based.',
        metavar='<number>'
    )
    p.add_argument(
        '--quant-cutoff',
        type=float,
        help='Quantification cutoff.',
        metavar='<number>'
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
    add_args_quiet(p)
    print_help_if_missing_args(p)
    p.set_defaults(func=filter_fasta)
    return p

def filter_fasta(args:argparse.Namespace) -> None:
    """ Filter noncanonical peptide FASTA """
    print_start_message(args)

    coding_tx = load_coding_transcripts(args)

    with open(args.input_fasta, 'rt') as handle:
        pool = VariantPeptidePool.load(handle)

    if not args.quiet:
        logger('Peptide FASTA file loaded.')

    with open(args.exprs_table, 'rt') as handle:
        i = 0
        while i < args.skip_lines:
            handle.readline()
            i += 1
        tx_id_col:str = args.tx_id_col
        quant_col:str = args.quant_col
        if any(not x.isdecimal() for x in [tx_id_col, quant_col]):
            header = handle.readline().rstrip().split(args.delimiter)

        if tx_id_col.isdecimal():
            tx_id_col = int(tx_id_col) - 1
        else:
            tx_id_col = header.index(tx_id_col)

        if quant_col.isdecimal():
            quant_col = int(quant_col) - 1
        else:
            quant_col = header.index(quant_col)

        exprs = load_expression_table(
            handle=handle,
            tx_col=tx_id_col,
            quant_col=quant_col,
            delim=args.delimiter
        )

    if not args.quiet:
        logger('Gene expression table loaded.')

    filtered_pool = pool.filter(exprs, args.quant_cutoff, coding_tx,
        args.keep_all_noncoding, args.keep_all_coding)

    filtered_pool.write(args.output_fasta)

    if not args.quiet:
        logger('Filtered FASTA file saved.')

def load_expression_table(handle:IO, tx_col:int,quant_col:int,
        delim:str='\t') -> Dict[str,float]:
    """ Load the gene expression quantification table """
    data = {}
    line:str
    for line in handle:
        fields = line.rstrip().split(delim)
        tx_id = fields[tx_col]
        quant = float(fields[quant_col])
        data[tx_id] = quant

    return data

def load_coding_transcripts(args:argparse.Namespace) -> List[str]:
    """ load and get the protein coding transcripts """
    if args.index_dir:
        with open(args.index_dir/'coding_transcripts.pkl', 'rb') as handle:
            return pickle.load(handle)

    anno = GenomicAnnotation()
    anno.dump_gtf(args.annotation_gtf)

    return [tx_id for tx_id, tx_model in anno.transcripts.items()
        if tx_model.is_protein_coding]
