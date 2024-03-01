""" Add fuzz test log """
import argparse
import datetime
import subprocess as sp
from typing import IO
from pathlib import Path
from moPepGen.cli.common import print_help_if_missing_args, add_args_debug_level
from moPepGen.util.fuzz_test import FuzzRecord, FuzzRecordStatus


FUZZ_TEST_LOG_HISTORY_CSV = Path(__file__).parent.parent.parent/'docs/files/fuzz_test_history.tsv'

def get_std(m:float, ss:float, n:int):
    """ Calculate standard deviation.

    Args:
        - m (float): mean
        - ss (float): sum of square
        - n (int): size
    """
    return (ss/n - m**2) ** (1/2)

# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ parse args """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='addFuzzTestLog',
        help='Add fuzz test log.'
    )
    p.add_argument(
        '-l', '--log-files',
        type=Path,
        help='Fuzz test log files',
        nargs='+',
        metavar='<files>',
        required=True
    )
    p.add_argument(
        '-c', '--commit-id',
        type=str,
        help='The Git commit ID that the fuzz tests used.',
        metavar='<value>',
        default=None
    )
    p.set_defaults(func=main)
    add_args_debug_level(p)
    print_help_if_missing_args(p)
    return p

class FuzzTestLogSummary:
    """ Fuzz test log summary """
    _commit = None
    _version = None

    def __init__(self, commit:str=None, version:str=None,
            n_match:int=0, n_mismatch:int=0, n_fail:int=0,
            avg_time_call_variant:datetime.timedelta=None,
            ss_time_call_variant:int=0,
            avg_time_brute_force:datetime.timedelta=None,
            ss_time_brute_force:int=0):
        """ constructor """
        self.commit = commit or self.get_commit()
        self.version = version or self.get_version()
        self.n_match = n_match
        self.n_mismatch = n_mismatch
        self.n_fail = n_fail
        self.avg_time_call_variant = avg_time_call_variant or datetime.timedelta()
        self.ss_time_call_variant = ss_time_call_variant
        self.avg_time_brute_force = avg_time_brute_force or datetime.timedelta()
        self.ss_time_brute_force = ss_time_brute_force

    @classmethod
    def get_commit(cls) -> str:
        """ Get the latest commit ID """
        if cls._commit:
            return cls._commit
        out = sp.check_output(['git', 'rev-parse', '--short', 'HEAD'])
        commit = out.decode('utf-8').rstrip()
        cls._commit = commit
        return commit

    @classmethod
    def get_version(cls) -> str:
        """ Get the latest version """
        if cls._version:
            return cls._version
        out = sp.check_output(['git', 'describe', '--tags', '--abbrev=0'])
        version = out.decode('utf-8').rstrip()
        cls._version = version
        return version

    def add_record(self, record:FuzzRecord):
        """ Add record """
        n_before = self.n_match + self.n_mismatch
        if record.status == FuzzRecordStatus.succeeded:
            self.n_match += 1
        elif record.status == FuzzRecordStatus.failed:
            self.n_mismatch += 1
            return
        else:
            self.n_fail += 1
            return

        dt_call_variant = (record.call_variant_end - record.call_variant_start)
        self.avg_time_call_variant = (
            self.avg_time_call_variant * n_before + dt_call_variant
        ) / (n_before + 1)
        self.ss_time_call_variant += dt_call_variant.total_seconds() ** 2

        dt_brute_force = (record.brute_force_end - record.brute_force_start)
        self.avg_time_brute_force = (
            self.avg_time_brute_force * n_before + dt_brute_force
        ) / (n_before + 1)
        self.ss_time_brute_force += dt_brute_force.total_seconds() ** 2


class FuzzTestLogSummaryCategorized:
    """ Fuzz test log summary separate by categories """
    def __init__(self, snv:FuzzTestLogSummary=None, indel:FuzzTestLogSummary=None,
            comprehensive:FuzzTestLogSummary=None, submit_date:datetime.date=None,
            commit:str=None):
        """ Constructor """
        self.snv = snv or FuzzTestLogSummary(commit=commit)
        self.indel = indel or FuzzTestLogSummary(commit=commit)
        self.comprehensive = comprehensive or FuzzTestLogSummary(commit=commit)
        self.submit_date = submit_date

    def summarize(self, handle:IO):
        """ Summarize fuzz run summary from the log file. """
        for record in FuzzRecord.parse(handle):
            if record.n_alt_splice == 0 \
                    and record.n_circ_rna == 0 \
                    and record.n_fusion == 0:
                if record.n_indel == 0:
                    self.snv.add_record(record)
                else:
                    self.indel.add_record(record)
            else:
                self.comprehensive.add_record(record)

            submitted = record.submitted.date()
            if self.submit_date is None or self.submit_date > submitted:
                self.submit_date = submitted

    @staticmethod
    def tsv_header() -> str:
        """ get TSV header """
        return '\t'.join([
            'version', 'commit', 'submit_date', 'variant_type',
            'n_match', 'n_mismatch', 'n_fail',
            'avg_time_call_variant', 'std_sec_call_variant',
            'avg_time_brute_force', 'std_sec_brute_force'
        ])

    def append_to_history(self, handle:IO):
        """ append summary to handle """
        if handle.tell() == 0:
            handle.write(self.tsv_header() + '\n')
        summaries = {
            'snv': self.snv,
            'indel': self.indel,
            'comprehensive': self.comprehensive
        }
        for key, summary in summaries.items():
            std_call_variant = get_std(
                m=summary.avg_time_call_variant.total_seconds(),
                ss=summary.ss_time_call_variant,
                n=summary.n_match
            )
            std_brute_force = get_std(
                m=summary.avg_time_brute_force.total_seconds(),
                ss=summary.ss_time_brute_force,
                n=summary.n_match
            )
            fields = [
                summary.version,
                summary.commit,
                str(self.submit_date),
                key,
                str(summary.n_match),
                str(summary.n_mismatch),
                str(summary.n_fail),
                str(summary.avg_time_call_variant),
                str(std_call_variant),
                str(summary.avg_time_brute_force),
                str(std_brute_force)
            ]
            handle.write('\t'.join(fields) + '\n')


def main(args:argparse.Namespace):
    """ Add fuzz test log """
    summary = FuzzTestLogSummaryCategorized(commit=args.commit_id)

    for file in args.log_files:
        with open(file, 'rt') as handle:
            summary.summarize(handle)

    with open(FUZZ_TEST_LOG_HISTORY_CSV, 'at') as handle:
        summary.append_to_history(handle)
