""" Module for fuzz test """
from __future__ import annotations
import argparse
import copy
from contextlib import redirect_stdout
from datetime import datetime, timedelta
from enum import Enum
import json
from pathlib import Path
import random
import shutil
import traceback
from typing import Iterable, List, Union, Tuple, IO
import sys
import uuid
from pathos.pools import ParallelPool
from Bio import SeqIO
from moPepGen import ERROR_INDEX_IN_INTRON, fake, seqvar, circ, dna, gtf, aa, get_logger
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.circ import CircRNAModel
from moPepGen.util.common import load_references
from moPepGen.cli.common import add_args_cleavage, generate_metadata, \
    print_help_if_missing_args, add_args_debug_level
from moPepGen.cli import call_variant_peptide
from moPepGen.gtf import GtfIO
from . import brute_force


class FuzzRecordStatus(str, Enum):
    """ Enumeration for fuzz record status. """
    succeeded = 'SUCCEEDED'
    failed = 'FAILED'

    def __str__(self) -> str:
        """ str """
        return self.value

# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction
        ) -> argparse.ArgumentParser:
    """ Arguments for the fuzzTest subcommand """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='fuzzTest',
        help='Fuzz test that randomly generates variants data and validate'
        ' the graph algorithm to brute force.'
    )
    parser.add_argument(
        '--reference-dir',
        type=Path,
        help='Directory to reference files. Noted this reference directory must'
        ' contain reference files generated by the downsampleReference command'
        ' and is different from the --index-dir for many moPepGen commands.',
        metavar = '<directory>',
        default=None
    )
    parser.add_argument(
        '--tx-id',
        type=str,
        help='Transcript ID',
        metavar='<value>'
    )
    parser.add_argument(
        '-s', '--max-size',
        type=int,
        help='Maximal number of nucleotides being inserted or deleted. Only'
        ' SNVs are generated when setting it to 0.',
        metavar='<number>'
    )
    parser.add_argument(
        '--max-variants',
        type=int,
        help='Max number of variants',
        metavar='<number>'
    )
    parser.add_argument(
        '--min-variants',
        type=int,
        help='Min number of variants',
        default=1,
        metavar='<number>'
    )
    parser.add_argument(
        '--snv-frac',
        type=float,
        default=1,
        help='Fraction of SNVs to simulate over INDEL in any test case. For'
        ' example, `--snv--frac 0.5` will have roughly 50% SNVs and 50% INDELs.'
        ' This has no effect on fusion, circRNA, or altSplice.'
    )
    parser.add_argument(
        '--snv-only-frac',
        type=float,
        default=1,
        help='Fraction of test cases to be SNV only (no INDEL). For example,'
        ' `--snv-only-frac 0.5` will have roughly 50% test cases to have only'
        ' SNV without any INDELs. This has no effect on fusion, circRNA, or altSplice.'
    )
    parser.add_argument(
        '--fusion-frac',
        type=float,
        help='Fraction of test cases to have fusion.',
        default=0
    )
    parser.add_argument(
        '--circ-rna-frac',
        type=float,
        help='Fraction of test cases to have circRNA.',
        default=0
    )
    parser.add_argument(
        '--ci-ratio',
        type=float,
        help='Fraction of ciRNA to generate.',
        default=0.2,
        metavar='<value>'
    )
    parser.add_argument(
        '--alt-splice-frac',
        type=float,
        help='Fraction of test cases to have alternative splicing.',
        default=0
    )
    parser.add_argument(
        '-i', '--n-iter',
        type=int,
        help='Max number of iterations',
        metavar='<number>'
    )
    parser.add_argument(
        '-t', '--temp-dir',
        type=Path,
        help='Temporary directory',
        metavar='<path>'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Output directory',
        metavar='<path>'
    )
    parser.add_argument(
        '--exonic-only',
        action='store_true',
        help='Create only exonic variants',
        default=False
    )
    parser.add_argument(
        '--keep-succeeded',
        action='store_true',
        help='Keep fuzz test files of the succeeded cases.'
    )
    parser.add_argument(
        '--save-graph',
        action='store_true',
        help='Save graph data.'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Random seed to set. If not specified, a random number will be generated.'
    )
    parser.add_argument(
        '--nthreads',
        type=int,
        default=1,
        help='Number of threads'
    )
    add_args_cleavage(parser)
    add_args_debug_level(parser)
    parser.set_defaults(func=main)
    print_help_if_missing_args(parser)
    return parser

def plain_logger(msg:str):
    """ log message to stdout """
    print(msg, flush=True, file=sys.stdout)

class FuzzRecord():
    """ Record of the fuzz test result for a single test case """
    def __init__(self, _id:str, status:str, submitted:datetime,
            completed:datetime, work_dir:Path=None, n_var:int=None,
            n_exonic:int=None, n_snv:int=0, n_indel:int=0, n_fusion:int=0,
            n_circ_rna:int=0, n_alt_splice:int=0, call_variant_start:datetime=None,
            call_variant_end:datetime=None, brute_force_start:datetime=None,
            brute_force_end:datetime=None):
        """ constructor """
        self.id = _id
        self.status = status
        self.submitted = submitted
        self.completed = completed
        self.work_dir = work_dir
        self.n_var = n_var
        self.n_exonic = n_exonic
        self.n_snv = n_snv
        self.n_indel = n_indel
        self.n_fusion = n_fusion
        self.n_circ_rna = n_circ_rna
        self.n_alt_splice = n_alt_splice
        self.n_snv = n_snv
        self.call_variant_start = call_variant_start
        self.call_variant_end = call_variant_end
        self.brute_force_start = brute_force_start
        self.brute_force_end = brute_force_end

    def submit(self):
        """ set the submit timestamp """
        self.submitted = datetime.now()

    def complete(self, status:str):
        """ set the complete timestamp """
        self.completed = datetime.now()
        self.status = status

    def call_variant_starts(self):
        """ set timestamp when callVariant starts """
        self.call_variant_start = datetime.now()

    def call_variant_ends(self):
        """ set timestamp when callVariant ends """
        self.call_variant_end = datetime.now()

    def brute_force_starts(self):
        """ set timestamp when bruteForce starts """
        self.brute_force_start = datetime.now()

    def brute_force_ends(self):
        """ set timestamp when bruteForce ends """
        self.brute_force_end = datetime.now()

    def set_work_dir(self, global_work_dir:Path):
        """ set test case specific work dir """
        self.work_dir = global_work_dir/self.id
        self.work_dir.mkdir(exist_ok=True)

    @classmethod
    def generate_record(cls) -> FuzzRecord:
        """ Generate a fuzz record """
        _id = uuid.uuid4().hex
        record =  cls(
            _id=_id, status=None, submitted=None, completed=None
        )
        record.submit()
        return record

    @classmethod
    def parse(cls, handle:IO) -> Iterable[FuzzRecord]:
        """ parse """
        handle.seek(0)
        header = handle.readline().rstrip().split('\t')
        time_fmt = '%Y-%m-%d %H:%M:%S'
        for line in handle:
            line:str
            fields = line.rstrip().split('\t')
            values = dict(zip(header, fields))

            if values['call_variant_start'] == '-':
                call_variant_start = None
            else:
                call_variant_start = datetime.strptime(values['call_variant_start'], time_fmt)

            if values['call_variant_end'] == '-':
                call_variant_end = None
            else:
                call_variant_end = datetime.strptime(values['call_variant_end'], time_fmt)

            if values['brute_force_start'] == '-':
                brute_force_start = None
            else:
                brute_force_start = datetime.strptime(values['brute_force_start'], time_fmt)

            if values['brute_force_end'] == '-':
                brute_force_end = None
            else:
                brute_force_end = datetime.strptime(values['brute_force_end'], time_fmt)

            yield FuzzRecord(
                _id=values['id'],
                status=values['status'],
                submitted=datetime.strptime(values['submitted'], time_fmt),
                completed=datetime.strptime(values['completed'], time_fmt),
                n_var=int(values['num_var']),
                n_exonic=int(values['num_exonic']),
                n_snv=int(values['n_snv']),
                n_indel=int(values['n_indel']),
                n_fusion=int(values['n_fusion']),
                n_circ_rna=int(values['n_circ_rna']),
                n_alt_splice=int(values['n_alt_splice']),
                call_variant_start=call_variant_start,
                call_variant_end=call_variant_end,
                brute_force_start=brute_force_start,
                brute_force_end=brute_force_end
            )

    @property
    def var_gvf_file(self) -> Path:
        """ File path to the GVF file that contains the faked variants """
        return self.work_dir/'fake_variants.gvf'

    @property
    def circ_gvf_file(self) -> Path:
        """ CircRNA GVF file """
        return self.work_dir/'fake_circ_rna.gvf'

    @property
    def call_variant_fasta(self) -> Path:
        """ FASTA file outputted by callVariant """
        return self.work_dir/'call_variant.fasta'

    @property
    def brute_force_fasta(self) -> Path:
        """ FASTA file outputted by bruteForce """
        return self.work_dir/'brute_force.fasta'

    @property
    def call_variant_only(self) -> Path:
        """ Path to the callVariant only file. """
        return self.work_dir/'call_variant_only.txt'

    @property
    def brute_force_only(self) -> Path:
        """ Path to the bruteForce only file. """
        return self.work_dir/'brute_force_only.txt'

    @property
    def timedelta(self) -> timedelta:
        """ Get time difference from submitted to completed """
        return self.completed - self.submitted

    def to_tsv(self) -> str:
        """ convert to TSV record """
        submitted = self.submitted.strftime('%Y-%m-%d %H:%M:%S') \
            if self.submitted else '-'
        completed = self.completed.strftime('%Y-%m-%d %H:%M:%S') \
            if self.completed else '-'
        call_variant_start = self.call_variant_start.strftime('%Y-%m-%d %H:%M:%S') \
            if self.call_variant_start else '-'
        call_variant_end = self.call_variant_end.strftime('%Y-%m-%d %H:%M:%S') \
            if self.call_variant_end else '-'
        call_variant_time = str(self.call_variant_end - self.call_variant_start) \
            if self.call_variant_start and self.call_variant_end else '-'
        brute_force_start = self.brute_force_start.strftime('%Y-%m-%d %H:%M:%S') \
            if self.brute_force_start else '-'
        brute_force_end = self.brute_force_end.strftime('%Y-%m-%d %H:%M:%S') \
            if self.brute_force_end else '-'
        brute_force_time = str(self.brute_force_end - self.brute_force_start) \
            if self.brute_force_start and self.brute_force_end else '-'
        return '\t'.join([
            self.id, self.status, str(self.n_var), str(self.n_exonic),
            str(self.n_snv), str(self.n_indel), str(self.n_fusion),
            str(self.n_circ_rna), str(self.n_alt_splice), submitted, completed,
            call_variant_start, call_variant_end, call_variant_time,
            brute_force_start, brute_force_end, brute_force_time
        ])

    @staticmethod
    def tsv_header() -> str:
        """ get TSV header """
        return '\t'.join([
            'id', 'status', 'num_var', 'num_exonic', 'n_snv', 'n_indel',
            'n_fusion', 'n_circ_rna', 'n_alt_splice', 'submitted', 'completed',
            'call_variant_start', 'call_variant_end', 'call_variant_time',
            'brute_force_start', 'brute_force_end', 'brute_force_time'
        ])

    @property
    def path_stderr(self) -> Path:
        """ File path to the test case specific stderr """
        return self.work_dir/'fuzzing.err'

    @property
    def path_stdout(self) -> Path:
        """ File path to the test case specific stdout """
        return self.work_dir/'fuzzing.out'

    def is_succeeded(self) -> bool:
        """ check if the test case was successful """
        return self.status == FuzzRecordStatus.succeeded

class FuzzTestCase():
    """ This models a single fuzz test case. It creates faked variants, save into
    a GVF file, and then calls moPepGen's callVariant and bruteForce in
    sequence. """
    def __init__(self, config:FuzzTestConfig, record:FuzzRecord=None):
        """ constructor """
        self.config = config
        self.record = record or FuzzRecord.generate_record()
        if not record:
            self.record.set_work_dir(self.config.temp_dir)

    def run(self):
        """ run the current test case """
        logger = get_logger()
        if self.config.seed is None:
            seed = random.randint(1, 1_000_000_000_000)
        else:
            seed = self.config.seed
        random.seed(seed)
        try:
            with open(self.record.path_stdout, 'w') as os_handle:
                with redirect_stdout(os_handle):
                    logger.info("[ fuzzTest ] Random seed: %i", seed)
                    self.fix_fusion_and_circ()
                    if self.config.ref_dir is None:
                        logger.info('[ fuzzTest ] Creating random reference.')
                        genome, anno, proteome = self.create_fake_reference()
                        self.write_fake_reference(genome, anno, proteome)
                        self.config.ref_dir = self.config.temp_dir/self.record.id
                        self.config.tx_id = list(anno.transcripts.keys())[0]
                    with open(self.record.work_dir/'config.json', 'w') as param_handle:
                        json.dump(self.config.to_dict(), param_handle, indent=4)
                    logger.info('Creating random variants.')
                    variants = self.create_variants()
                    self.write_variants(variants)
                    self.record.call_variant_starts()
                    logger.info('Starting callVariant.')
                    self.call_variants()
                    self.record.call_variant_ends()
                    self.record.brute_force_starts()
                    logger.info('Starting bruteForce.')
                    self.brute_force()
                    self.record.brute_force_ends()
                    status = self.assert_equal()
                    self.record.complete(status)
            if status == FuzzRecordStatus.succeeded and not self.config.keep_succeeded:
                shutil.rmtree(self.config.temp_dir/self.record.id)
            else:
                shutil.copytree(
                    self.config.temp_dir/self.record.id,
                    self.config.output_dir/self.record.id
                )
        # pylint: disable=W0703
        except Exception as error:
            self.record.complete(str(FuzzRecordStatus.failed))
            with self.record.path_stderr.open('wt') as handle:
                handle.write(str(error))
                handle.write(traceback.format_exc())
            shutil.copytree(
                self.config.temp_dir/self.record.id,
                self.config.output_dir/self.record.id
            )

    def fix_fusion_and_circ(self):
        """ Fix number of fusion and circRNA """
        if random.choices(
            [True, False],
            weights=[self.config.fusion_frac, 1-self.config.fusion_frac]
        )[0]:
            self.record.n_fusion = 1

        if random.choices(
            [True, False],
            weights=[self.config.circ_rna_frac, 1-self.config.circ_rna_frac]
        )[0]:
            self.record.n_circ_rna = 1

    def create_fake_reference(self
            ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation, aa.AminoAcidSeqDict]:
        """ Create a fake genomic reference, including the genome sequence and
        genomic annotation model, and the translated amino acid sequences. """
        if self.record.n_fusion > 0:
            n_genes = 2
        else:
            n_genes = 1
        genome, anno = fake.fake_genome_and_annotation(n_genes)
        proteome = aa.AminoAcidSeqDict()

        for tx_model in anno.transcripts.values():
            if not tx_model.is_protein_coding:
                continue
            tx_seq = tx_model.get_transcript_sequence(genome[tx_model.transcript.chrom])
            aa_seq = tx_seq[tx_seq.orf.start:tx_seq.orf.end].translate()
            for sec in tx_seq.selenocysteine:
                sec_start = int((sec.start - tx_seq.orf.start) / 3)
                aa_seq = aa_seq[:sec_start] + 'U' + aa_seq[sec_start+1:]
            aa_seq.description = \
                f"{tx_model.protein_id}|{tx_model.transcript_id}|{tx_model.gene_id}|XXX"
            proteome[tx_model.transcript_id] = aa_seq

        return genome, anno, proteome

    def write_fake_reference(self, genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation,
            proteome:aa.AminoAcidSeqDict):
        """ Write the fake genomic reference into files. """
        work_dir = self.config.temp_dir/self.record.id

        with open(work_dir/'annotation.gtf', 'wt') as handle:
            GtfIO.write(handle, anno)

        with open(work_dir/'genome.fasta', 'wt') as handle:
            writer = SeqIO.FastaIO.FastaWriter(handle, record2title=lambda x:x.id)
            for record in genome.values():
                writer.write_record(record)

        with open(work_dir/'proteome.fasta', 'wt') as handle:
            writer = SeqIO.FastaIO.FastaWriter(handle, record2title=lambda x: x.description)
            for record in proteome.values():
                writer.write_record(record)

    def create_variants(self) -> List[Union[VariantRecord,CircRNAModel]]:
        """ create random variant records """
        anno, genome, _ = load_references(
            path_anno=self.config.path_annotation_gtf,
            path_genome=self.config.path_genome_fasta,
            path_proteome=self.config.path_proteome_fasta
        )
        records:List[Union[VariantRecord, CircRNAModel]] = []
        n_variants = random.randint(self.config.min_vairants, self.config.max_variants)
        tx_ids = [self.config.tx_id]

        if self.record.n_fusion > 0:
            if random.choice([True, False]):
                self.record.n_fusion = 1
                record = fake.fake_fusion(anno, genome, self.config.tx_id)
                records.append(record)
                n_variants -= 1
                tx_ids.append(record.attrs['ACCEPTER_TRANSCRIPT_ID'])

        if self.record.n_circ_rna > 0:
            if random.choice([True, False]):
                self.record.n_circ_rna = 1
                record = fake.fake_circ_rna_model(anno, self.config.tx_id,
                    self.config.ci_ratio)
                records.append(record)
                n_variants -= 1

        var_types = ['SNV']
        if random.choices(
            [True, False],
            weights=[self.config.alt_splice_frac, 1-self.config.alt_splice_frac]
        )[0]:
            var_types.append('AltSplicing')
        snv_only = random.choices(
            [True, False],
            weights=[self.config.snv_only_frac, 1-self.config.snv_only_frac]
        )[0]
        for _ in range(n_variants):
            tx_id = random.choice(tx_ids)
            var_type = random.choice(var_types)
            if var_type == 'AltSplicing':
                record = fake.fake_rmats_record(anno, genome, tx_id)
                self.record.n_alt_splice += 1
            elif var_type == 'SNV':
                if snv_only:
                    var_type_sub = 'SNV'
                else:
                    var_type_sub = random.choices(
                        ['SNV', 'INDEL'],
                        weights=[self.config.snv_frac, 1-self.config.snv_frac]
                    )[0]
                record = fake.fake_variant_record(
                    anno=anno, genome=genome, tx_id=tx_id,
                    var_type=var_type_sub,
                    max_size=self.config.max_size,
                    exonic_only=self.config.exonic_only
                )
                if record.type == 'SNV':
                    self.record.n_snv += 1
                else:
                    self.record.n_indel += 1
            else:
                raise ValueError(f"var_type of {var_type} can not be recognized.")
            records.append(record)
        self.record.n_var = n_variants
        self.record.n_exonic = 0
        for variant in records:
            if variant.__class__ is CircRNAModel:
                continue
            try:
                variant.to_transcript_variant(anno, genome, variant.transcript_id)
                self.record.n_exonic += 1
            except ValueError as e:
                if e.args[0] == ERROR_INDEX_IN_INTRON:
                    continue
                raise e

        return records

    def write_variants(self, variants:Iterable[Union[VariantRecord,CircRNAModel]]):
        """ Write variants to file """
        var_records = []
        circ_records = []

        for v in variants:
            if v.__class__ is CircRNAModel:
                circ_records.append(v)
            else:
                var_records.append(v)

        args = argparse.Namespace()
        args.index_dir = None
        args.command = 'parseVEP'
        args.source = ''
        metadata = generate_metadata(args)
        seqvar.io.write(var_records, self.record.var_gvf_file, metadata)

        args.command = 'parseCIRCexplorer'
        metadata = generate_metadata(args)
        if circ_records:
            with open(self.record.circ_gvf_file, 'w') as handle:
                circ.io.write(circ_records, metadata, handle)

    def call_variants(self):
        """ call variants using moPepGen's graph algorithm """
        args = argparse.Namespace()
        args.index_dir = None
        args.command = 'callPeptides'
        args.input_path = [self.record.var_gvf_file]
        if self.record.circ_gvf_file.exists():
            args.input_path.append(self.record.circ_gvf_file)
        args.genome_fasta = self.config.path_genome_fasta
        args.annotation_gtf = self.config.path_annotation_gtf
        args.proteome_fasta = self.config.path_proteome_fasta
        args.reference_source = None
        args.output_path = self.record.call_variant_fasta
        if self.config.save_graph:
            args.graph_output_dir = self.record.work_dir/'graph'
            args.graph_output_dir.mkdir(exist_ok=True)
        else:
            args.graph_output_dir = None
        args.quiet = True
        args.backsplicing_only = False
        args.max_adjacent_as_mnv = 2
        args.selenocysteine_termination = True
        args.w2f_reassignment = True
        args.cleavage_rule = 'trypsin'
        args.cleavage_exception = None
        args.miscleavage = 2
        args.min_mw = 500.
        args.min_length = 7
        args.max_length = 25
        args.threads = 1
        args.max_variants_per_node = (9, )
        args.additional_variants_per_misc = (2, )
        args.min_nodes_to_collapse = 30
        args.naa_to_collapse = 5
        args.noncanonical_transcripts = False
        args.invalid_protein_as_noncoding = False
        args.debug_level = 1
        args.timeout_seconds = 300
        args.coding_novel_orf = False
        args.skip_failed = False
        call_variant_peptide(args=args)

    def brute_force(self):
        """ call the brute force variant peptide caller """
        args = argparse.Namespace()
        args.input_gvf = [self.record.var_gvf_file]
        if self.record.circ_gvf_file.exists():
            args.input_gvf.append(self.record.circ_gvf_file)
        args.reference_dir = self.config.ref_dir
        args.force = True
        args.variant_ids = []
        args.max_adjacent_as_mnv = 2
        args.selenocysteine_termination = True
        args.w2f_reassignment = True
        args.cleavage_rule = self.config.cleavage_rule
        args.cleavage_exception = None
        args.miscleavage = self.config.miscleavage
        args.min_mw = self.config.min_mw
        args.min_length = self.config.min_length
        args.max_length = self.config.max_length

        with open(self.record.brute_force_fasta, 'wt') as handle:
            with redirect_stdout(handle):
                brute_force.main(args)

    def assert_equal(self) -> str:
        """ Assert that the callVariant results and bruteForce results equal """
        with open(self.record.call_variant_fasta, 'rt') as handle:
            variant_seqs = set()
            for seq in SeqIO.parse(handle, 'fasta'):
                variant_seqs.add(str(seq.seq))
        with open(self.record.brute_force_fasta, 'rt') as handle:
            brute_force_seqs = set()
            for line in handle:
                brute_force_seqs.add(line.rstrip())

        if variant_seqs != brute_force_seqs:
            variant_only = variant_seqs - brute_force_seqs
            brute_force_only = brute_force_seqs - variant_seqs
            if variant_only:
                with open(self.record.call_variant_only, 'wt') as handle:
                    for seq in variant_only:
                        handle.write(seq + '\n')
            if brute_force_only:
                with open(self.record.brute_force_only, 'wt') as handle:
                    for seq in brute_force_only:
                        handle.write(seq + '\n')
            return str(FuzzRecordStatus.failed)
        return str(FuzzRecordStatus.succeeded)

class FuzzTestConfig():
    """ Fuzz test config """
    def __init__(self, tx_id:str, n_iter:int, max_size:int, max_variants:int,
            min_variants:int, exonic_only:bool, snv_frac:float, snv_only_frac:float,
            fusion_frac:float, circ_rna_frac:float, ci_ratio:float, alt_splice_frac:bool,
            keep_succeeded:bool, save_graph:bool, cleavage_rule:str, cleavage_exception:str,
            miscleavage:int, min_mw:int, min_length:int, max_length:int, temp_dir:Path,
            output_dir:Path, ref_dir:Path, fuzz_start:datetime=None, fuzz_end:datetime=None,
            seed:int=None, nthreads:int=1):
        """ constructor """
        self.tx_id = tx_id
        self.n_iter = n_iter
        self.max_size = max_size
        self.max_variants = max_variants
        self.min_vairants = min_variants
        self.exonic_only = exonic_only
        self.snv_frac = snv_frac
        self.snv_only_frac = snv_only_frac
        self.fusion_frac = fusion_frac
        self.circ_rna_frac = circ_rna_frac
        self.ci_ratio = ci_ratio
        self.alt_splice_frac = alt_splice_frac
        self.keep_succeeded = keep_succeeded
        self.save_graph = save_graph
        self.cleavage_rule = cleavage_rule
        self.cleavage_exception = cleavage_exception
        self.miscleavage = miscleavage
        self.min_mw = min_mw
        self.min_length = min_length
        self.max_length = max_length
        self.temp_dir = temp_dir
        self.output_dir = output_dir
        self.ref_dir = ref_dir
        self.fuzz_start = fuzz_start
        self.fuzz_end = fuzz_end
        self.seed = seed
        self.nthreads = nthreads

    def to_dict(self):
        """ Convert to a dict """
        mapper = {}
        for k,v in self.__dict__.items():
            if isinstance(v, Path):
                v = str(v)
            elif isinstance(v, datetime):
                continue
            mapper[k] = v
        return mapper

    @property
    def path_genome_fasta(self) -> Path:
        """ Genome FASTA path """
        return self.ref_dir/'genome.fasta'

    @property
    def path_annotation_gtf(self) -> Path:
        """ Annotation FASTA path """
        return self.ref_dir/'annotation.gtf'

    @property
    def path_proteome_fasta(self) -> Path:
        """ Proteome FASTA path """
        return self.ref_dir/'proteome.fasta'

    @property
    def path_log_file(self) -> Path:
        """ Log file path """
        start_time = self.fuzz_start.strftime('%y%m%d-%H%M%S')
        return self.output_dir/f"fuzz-{start_time}.log"

class Fuzzer():
    """ Main class for fuzzing """
    def __init__(self, config:FuzzTestConfig, records:List[FuzzRecord]=None):
        """ constructor """
        self.config = config
        self.config.temp_dir.mkdir(exist_ok=True)
        self.config.output_dir.mkdir(exist_ok=True)
        self.records = records or []
        self.handle = None

    def __enter__(self):
        """ enter """
        self.set_fuzz_start_time()
        self.print_prolog()
        self.handle = self.config.path_log_file.open('wt')
        self.handle.write(FuzzRecord.tsv_header() + '\n')
        self.handle.flush()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ exit """
        self.handle.close()
        self.set_fuzz_end_time()
        self.print_epilog()

    def set_fuzz_start_time(self):
        """ set fuzz start time """
        self.config.fuzz_start = datetime.now()

    def set_fuzz_end_time(self):
        """ set fuzz end time """
        self.config.fuzz_end = datetime.now()

    def print_prolog(self):
        """ print prolog """
        fuzz_start = self.config.fuzz_start.strftime('%Y-%m-%d %H:%M:%S')
        plain_logger(f"moPepGen fuzz testing started at {fuzz_start}")

    def print_epilog(self):
        """ print epilog """
        fuzz_end = self.config.fuzz_end.strftime('%Y-%m-%d %H:%M:%S')
        plain_logger(f"moPepGen fuzz testing ended at {fuzz_end}")
        n_succeeded = 0
        n_failed = 0
        for x in self.records:
            if x.is_succeeded():
                n_succeeded += 1
            else:
                n_failed += 1
        plain_logger(f"{n_succeeded} cases passed.")
        plain_logger(f"{n_failed} cases failed.")

    def fuzz(self):
        """ fuzz """
        pool = ParallelPool(ncpus=self.config.nthreads)
        cases:List[FuzzTestCase] = []
        for i in range(self.config.n_iter):
            config = copy.deepcopy(self.config)
            case = FuzzTestCase(config)
            cases.append(case)
            if len(cases) >= config.nthreads or i == self.config.n_iter - 1:
                res:List[FuzzTestCase]
                if self.config.nthreads == 1:
                    res = map(_run_fuzz, cases)
                else:
                    res = pool.map(_run_fuzz, cases)
                for case in res:
                    self.handle.write(case.record.to_tsv() + '\n')
                    self.handle.flush()
                    self.records.append(case.record)
                cases = []

def _run_fuzz(case:FuzzTestCase) -> FuzzTestCase:
    """ Run fuzz test """
    case.run()
    return case

def main(args:argparse.Namespace):
    """ Main entry point for fuzz test """
    config = FuzzTestConfig(
        tx_id=args.tx_id,
        n_iter=args.n_iter,
        max_size=args.max_size,
        max_variants=args.max_variants,
        min_variants=args.min_variants,
        exonic_only=args.exonic_only,
        snv_frac=args.snv_frac,
        snv_only_frac=args.snv_only_frac,
        fusion_frac=args.fusion_frac,
        alt_splice_frac=args.alt_splice_frac,
        circ_rna_frac=args.circ_rna_frac,
        ci_ratio=args.ci_ratio,
        cleavage_rule=args.cleavage_rule,
        cleavage_exception=args.cleavage_exception,
        miscleavage=args.miscleavage,
        min_mw=args.min_mw,
        min_length=args.min_length,
        max_length=args.max_length,
        temp_dir=args.temp_dir,
        output_dir=args.output_dir,
        ref_dir=args.reference_dir,
        keep_succeeded=args.keep_succeeded,
        save_graph=args.save_graph,
        seed=args.seed,
        nthreads=args.nthreads
    )
    with Fuzzer(config) as fuzzer:
        fuzzer.fuzz()
