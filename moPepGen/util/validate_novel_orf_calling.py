""" Validate callNovelORF results with bruteForceNoncoding """
import argparse
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
import random
from typing import IO, List, Set
from Bio import SeqIO
from moPepGen import get_logger
from moPepGen.cli import call_novel_orf_peptide, common
from moPepGen.gtf import GtfIO
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen.util.common import load_references
from moPepGen.util.validate_variant_calling import call_downsample_reference
from moPepGen.util import brute_force_novel_orf


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='validateNovelORFCalling',
        help='Validate the novel ORF peptide calling result of the graph-based'
        ' algorithm with the brute force algorithm.'
    )
    parser.add_argument(
        '--tx-id',
        type=str,
        nargs='*',
        metavar='<value>',
        help='Transcript IDs to sample from to validate. If not given, all'
        ' noncoding transcripts from the annotation GTF will be used.'
    )
    parser.add_argument(
        '--n-iter',
        type=int,
        help='Number of iterations',
        metavar='<number>',
        default=1
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Path to the output diractory.',
        metavar='<file>'
    )
    common.add_args_reference(parser, proteome=True, index=False)
    parser.set_defaults(func=main)
    common.print_help_if_missing_args(parser)
    return parser

def get_transcript_ids(handle:IO) -> Set[str]:
    """ Get all transcript IDs. """
    _, exclusion_biotypes = common.load_inclusion_exclusion_biotypes(
        argparse.Namespace(inclusion_biotypes=None, exclusion_biotypes=None)
    )
    tx_ids:Set[str]= set()
    record:GTFSeqFeature
    for record in GtfIO.parse(handle):
        record.source = 'GENCODE'
        if record.type.lower() == 'gene':
            continue
        if record.biotype in exclusion_biotypes:
            continue
        tx_ids.add(record.transcript_id)
    return tx_ids

def call_novel_orf(ref_dir:Path, output_fasta:Path):
    """ call callNovelORF """
    args = argparse.Namespace()
    args.command = 'callNovelORF'
    args.genome_fasta = ref_dir/'genome.fasta'
    args.annotation_gtf = ref_dir/'annotation.gtf'
    args.proteome_fasta = ref_dir/'proteome.fasta'
    args.index_dir = None
    args.min_tx_length = 21
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25
    args.output_orf = None
    args.output_path = output_fasta
    args.min_tx_length = 21
    args.coding_novel_orf = False
    args.inclusion_biotypes = None
    args.exclusion_biotypes = None
    args.quiet = True
    call_novel_orf_peptide(args)

def call_brute_force_noncoding(tx_id:str, ref_dir:Path, output_path:Path):
    """ call bruteForceNoncoding """
    args = argparse.Namespace()
    args.tx_id = tx_id
    args.reference_dir = ref_dir
    args.canonical_peptides = None
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    args.min_length = 7
    args.max_length = 25

    with open(output_path, 'wt') as handle:
        with redirect_stdout(handle):
            brute_force_novel_orf.main(args)

def get_transcript_length(tx_id:str, ref_dir:Path) -> int:
    """ Get the transcript length """
    anno, genome, *_ = load_references(
        path_anno=ref_dir/'annotation.gtf',
        path_genome=ref_dir/'genome.fasta',
        path_proteome=ref_dir/'proteome.fasta'
    )
    return len(anno.transcripts[tx_id].get_transcript_sequence(genome['chr1']))

def assert_equal(noncoding_fasta:Path, brute_force_txt:Path, output_dir:Path) -> bool:
    """ Compare the results of callNovelORF and bruteForceNovelORF """
    with open(noncoding_fasta, 'rt') as handle:
        noncoding_seqs = set()
        for seq in SeqIO.parse(handle, 'fasta'):
            noncoding_seqs.add(str(seq.seq))

    with open(brute_force_txt, 'rt') as handle:
        brute_force_seqs = set()
        for line in handle:
            brute_force_seqs.add(line.rstrip())

    if noncoding_seqs != brute_force_seqs:
        noncoding_only = noncoding_seqs - brute_force_seqs
        brute_force_only = brute_force_seqs - noncoding_seqs

        if noncoding_only:
            noncoding_only_file = output_dir/'call_noncoding_only.txt'
            with open(noncoding_only_file, 'wt') as handle:
                for seq in noncoding_only:
                    handle.write(seq + '\n')

        if brute_force_only:
            brute_force_only_file = output_dir/'brute_force_only.txt'
            with open(brute_force_only_file, 'wt') as handle:
                for seq in brute_force_only:
                    handle.write(seq + '\n')

        return False
    return True

class ValidationRecord():
    """ Validation summary record """
    def __init__(self, tx_id:str=None, status:str=None, submitted:datetime=None,
            completed:datetime=None, tx_len:int=None,
            call_noncoding_start:datetime=None, call_noncoding_end:datetime=None,
            brute_force_start:datetime=None, brute_force_end:datetime=None):
        """ """
        self.tx_id = tx_id
        self.status = status
        self.submitted = submitted
        self.completed = completed
        self.tx_len = tx_len
        self.call_noncoding_start = call_noncoding_start
        self.call_noncoding_end = call_noncoding_end
        self.brute_force_start = brute_force_start
        self.brute_force_end = brute_force_end

    def submit(self):
        """ Validation starts """
        self.submitted = datetime.now()

    def complete(self, status:str):
        """ set the complete timestamp """
        self.completed = datetime.now()
        self.status = status

    def start_call_novel_orf(self):
        """ set timestamp when callNovelORF starts """
        self.call_noncoding_start = datetime.now()

    def end_call_novel_orf(self):
        """ set timestamp when callNovelORF ends """
        self.call_noncoding_end = datetime.now()

    def start_brute_force(self):
        """ set timestamp when bruteForce starts """
        self.brute_force_start = datetime.now()

    def end_brute_force(self):
        """ set timestamp when bruteForce ends """
        self.brute_force_end = datetime.now()

    def to_tsv(self) -> str:
        """ convert to TSV record """
        submitted = self.submitted.strftime('%Y-%m-%d %H:%M:%S') \
            if self.submitted else '-'
        completed = self.completed.strftime('%Y-%m-%d %H:%M:%S') \
            if self.completed else '-'
        call_noncoding_start = self.call_noncoding_start.strftime('%Y-%m-%d %H:%M:%S') \
            if self.call_noncoding_start else '-'
        call_noncoding_end = self.call_noncoding_end.strftime('%Y-%m-%d %H:%M:%S') \
            if self.call_noncoding_end else '-'
        call_noncoding_time = str(self.call_noncoding_end - self.call_noncoding_start) \
            if self.call_noncoding_start and self.call_noncoding_end else '-'
        brute_force_start = self.brute_force_start.strftime('%Y-%m-%d %H:%M:%S') \
            if self.brute_force_start else '-'
        brute_force_end = self.brute_force_end.strftime('%Y-%m-%d %H:%M:%S') \
            if self.brute_force_end else '-'
        brute_force_time = str(self.brute_force_end - self.brute_force_start) \
            if self.brute_force_start and self.brute_force_end else '-'
        return '\t'.join([
            self.tx_id, str(self.tx_len), self.status, submitted, completed,
            call_noncoding_start, call_noncoding_end, call_noncoding_time,
            brute_force_start, brute_force_end, brute_force_time
        ])


class ValidationSummary():
    """ Valiadtion summary """
    def __init__(self, data:List[ValidationRecord]=None):
        """ constructor """
        self.data = data or []

    def append(self, x:ValidationRecord):
        """ Add summary record """
        self.data.append(x)

    def write(self, path:Path):
        """ Write summary """
        header = '\t'.join([
            'tx_id', 'tx_len', 'status',
            'submitted', 'completed',
            'call_noncoding_start', 'call_noncoding_end', 'call_noncoding_time',
            'brute_force_start', 'brute_force_end', 'brute_force_time'
        ])
        with open(path, 'w') as handle:
            handle.write(header + '\n')
            for record in self.data:
                handle.write(record.to_tsv() + '\n')

def main(args:argparse.Namespace):
    """ main entrypoint """
    logger = get_logger()
    if args.tx_id:
        tx_ids = args.tx_id
    else:
        with open(args.annotation_gtf, 'rt') as handle:
            tx_ids = get_transcript_ids(handle)

    n_iter = min(args.n_iter, len(tx_ids))
    tx_ids_sampled = random.sample(tx_ids, n_iter)

    summary = ValidationSummary()

    for tx_id in tx_ids_sampled:
        output_dir:Path = args.output_dir/tx_id
        output_dir.mkdir(exist_ok=True)
        ref_dir = output_dir/'index'
        ref_dir.mkdir(exist_ok=True)

        record = ValidationRecord(tx_id=tx_id)
        summary.append(record)
        record.submit()

        noncoding_fasta = output_dir/'call_noncoding.fasta'
        brute_force_txt = output_dir/'brute_force_noncoding.fasta'

        call_downsample_reference(
            genome=args.genome_fasta,
            anno=args.annotation_gtf,
            protein=args.proteome_fasta,
            tx_id=tx_id,
            output_dir=ref_dir
        )

        record.tx_len = get_transcript_length(tx_id, ref_dir)

        record.start_call_novel_orf()
        call_novel_orf(
            ref_dir=ref_dir,
            output_fasta=noncoding_fasta
        )
        record.end_call_novel_orf()

        record.start_brute_force()
        call_brute_force_noncoding(
            tx_id=tx_id,
            ref_dir=ref_dir,
            output_path=brute_force_txt
        )
        record.end_brute_force()

        res = assert_equal(
            noncoding_fasta=noncoding_fasta,
            brute_force_txt=brute_force_txt,
            output_dir=output_dir
        )
        record.complete('SUCCEEDED' if res else 'FAILED')

        logger.info(
            "Transcript ID: %s, %s!",
            tx_id,
            'Equal' if res else 'Not equal'
        )

    summary.write(args.output_dir/'validate_noncoding_summary.tsv')
