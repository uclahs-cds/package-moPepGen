""" `filterFasta` takes the FASTA file of variant peptides (output by
`callVariant`) or novel ORF peptides (output by `callNovelORF`) and filters it
based on the gene expression data. A expresion table must be given as a CSV
or TSV. """
from __future__ import annotations
import argparse
from pathlib import Path
import pickle
from typing import IO, Dict, List
from Bio import SeqIO
from moPepGen import get_logger
from moPepGen.aa import VariantPeptidePool
from moPepGen.aa.expasy_rules import EXPASY_RULES
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation
from moPepGen.cli import common


INPUT_FILE_FORMATS = ['.fasta', '.fa']
OUTPUT_FILE_FORMATS = ['.fasta', '.fa']

# pylint: disable=W0212
def add_subparser_filter_fasta(subparser:argparse._SubParsersAction):
    """ CLI for moPepGen filterFasta """
    p:argparse.ArgumentParser = subparser.add_parser(
        name='filterFasta',
        help='Filter noncanonical peptides.',
        description='Filter noncanonical peptides according to gene expression'
        ' or gene biotypes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=p, formats=INPUT_FILE_FORMATS,
        message='Input FASTA file, must be generated by either moPepGen'
        ' callVariant or callNovelORF.'
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--denylist',
        help='Path to the peptide sequence deny list. When using novel ORF'
        ' peptides as denylist, make sure it is no also passed as a input FASTA'
        ' file, because all peptide sequences will be removed. Valid formats:'
        ' [".fasta"]',
        type=Path,
        metavar='<file>',
        default=None
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
        help='The index for transcript ID in the RNAseq quantification results.'
        ' Index is 1-based.',
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
        metavar='<number>',
        default=None
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
    p.add_argument(
        '--keep-canonical',
        action='store_true',
        help='Keep peptides called from canonical ORFs. Only useful together'
        ' with denylist.'
    )
    p.add_argument(
        '--miscleavages',
        type=str,
        help='Range of miscleavages per peptide to allow. Min and max are'
        ' included. For example, "1:2" will keep all peptides with 1 or 2'
        ' miscleavages.',
        metavar='[min]:[max]',
        default=None
    )
    p.add_argument(
        '--enzyme',
        type=str,
        help='Enzyme name. Ignored if --miscleavages is not specified.',
        choices=list(EXPASY_RULES.keys()),
        default='trypsin'
    )

    common.add_args_reference(p, genome=False, proteome=False)
    common.add_args_debug_level(p)
    common.print_help_if_missing_args(p)
    p.set_defaults(func=filter_fasta)
    return p

def filter_fasta(args:argparse.Namespace) -> None:
    """ Filter noncanonical peptide FASTA """
    logger = get_logger()
    common.validate_file_format(
        args.input_path, INPUT_FILE_FORMATS, check_readable=True
    )
    common.validate_file_format(
        args.output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    if args.miscleavages is not None and not ':' in args.miscleavages:
        raise ValueError('Value invalid for --miscleavages. Must be [min]:[max]')

    if args.miscleavages is None:
        miscleavage_range = (None, None)
    else:
        miscleavage_range = tuple(int(x) for x in args.miscleavages.split(':', 1))

    common.print_start_message(args)

    coding_tx = load_coding_transcripts(args)

    with open(args.input_path, 'rt') as handle:
        pool = VariantPeptidePool.load(handle)

    logger.info('Peptide FASTA file loaded.')

    if args.exprs_table is not None:
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
        logger.info('Gene expression table loaded.')
    else:
        exprs = None

    if args.denylist is not None:
        with open(args.denylist, 'rt') as handle:
            denylist = {seq.seq for seq in SeqIO.parse(handle, 'fasta')}
        logger.info('Peptide denylist loaded.')
    else:
        denylist = None

    filtered_pool = pool.filter(
        exprs=exprs, cutoff=args.quant_cutoff, coding_transcripts=coding_tx,
        keep_all_noncoding=args.keep_all_noncoding,
        keep_all_coding=args.keep_all_coding, enzyme=args.enzyme,
        miscleavage_range=miscleavage_range, denylist=denylist,
        keep_canonical=args.keep_canonical
    )

    filtered_pool.write(args.output_path)
    logger.info('Filtered FASTA file saved.')

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

    return {tx_id for tx_id, tx_model in anno.transcripts.items()
        if tx_model.is_protein_coding}
