""" Add fuzz test log """
import argparse
import datetime
import subprocess as sp
from typing import IO
from pathlib import Path
from moPepGen.cli.common import print_help_if_missing_args
from moPepGen.util.fuzz_test import FuzzRecord, FuzzRecordStatus


FUZZ_TEST_LOG_HISTORY_CSV = Path(__file__).parent.parent.parent/'docs/files/fuzz_test_history.tsv'

# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
    """ """
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
    print_help_if_missing_args(p)
    return p

class FuzzTestLogSummary:
    """ Fuzz test log summary """
    _commit = None
    _version = None

    def __init__(self, commit:str=None, version:str=None,
            n_match:int=0, n_mismatch:int=0, n_fail:int=0,
            avg_time_call_variant:datetime.timedelta=None,
            avg_time_brute_force:datetime.timedelta=None):
        """ constructor """
        self.commit = commit or self.get_commit()
        self.version = version or self.get_version()
        self.n_match = n_match
        self.n_mismatch = n_mismatch
        self.n_fail = n_fail
        self.avg_time_call_variant = avg_time_call_variant or datetime.timedelta()
        self.avg_time_brute_force = avg_time_brute_force or datetime.timedelta()

    @classmethod
    def get_commit(cls) -> str:
        """ """
        if cls._commit:
            return cls._commit
        out = sp.check_output(['git', 'rev-parse', '--short', 'HEAD'])
        commit = out.decode('utf-8').rstrip()
        cls._commit = commit
        return commit

    @classmethod
    def get_version(cls) -> str:
        """ """
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
        else:
            self.n_fail += 1
            return

        self.avg_time_call_variant = (
            self.avg_time_call_variant * n_before
            + (record.call_variant_end - record.call_variant_start)
        ) / (n_before + 1)
        self.avg_time_brute_force = (
            self.avg_time_brute_force * n_before
            + (record.brute_force_end - record.brute_force_start)
        ) / (n_before + 1)


class FuzzTestLogSummaryCategorized:
    """ Fuzz test log summary separate by categories """
    def __init__(self, snv:FuzzTestLogSummary=None, indel:FuzzTestLogSummary=None,
            comprehensive:FuzzTestLogSummary=None, submit_date:datetime.date=None,
            commit:str=None):
        """ """
        self.snv = snv or FuzzTestLogSummary(commit=commit)
        self.indel = indel or FuzzTestLogSummary(commit=commit)
        self.comprehensive = comprehensive or FuzzTestLogSummary(commit=commit)
        self.submit_date = submit_date

    def summarize(self, handle:IO):
        """ """
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
            'avg_time_call_variant', 'avg_time_brute_force'
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
            fields = [
                summary.version,
                summary.commit,
                str(self.submit_date),
                key,
                str(summary.n_match),
                str(summary.n_mismatch),
                str(summary.n_fail),
                str(summary.avg_time_call_variant),
                str(summary.avg_time_brute_force)
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
