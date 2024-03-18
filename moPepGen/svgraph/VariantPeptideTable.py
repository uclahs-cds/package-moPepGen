""" Varaint Peptide Table """
from __future__ import annotations
from typing import TYPE_CHECKING, Dict, Set, IO, List, Tuple
from Bio import SeqUtils
from Bio.SeqIO import FastaIO
from moPepGen import VARIANT_PEPTIDE_SOURCE_DELIMITER, aa
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.svgraph.VariantPeptideDict import AnnotatedPeptideLabel, PeptideSegment


if TYPE_CHECKING:
    from pathlib import Path
    from Bio.Seq import Seq
    from moPepGen.params import CleavageParams

VARIANT_PEPTIDE_TABLE_HEADERS = [
    'sequence', 'header', 'subsequence', 'start', 'end', 'feature_type', 'feature_id',
    'ref_start', 'ref_end', 'start_offset', 'end_offset', 'variant'
]

class VariantPeptideTable:
    """ Variant Peptide Segment Annotation """
    def __init__(self, handle:IO[str], index:Dict[Seq, List[Tuple[int,int]]]=None):
        self.handle = handle
        self.index = index or {}
        self.header_delimeter = VARIANT_PEPTIDE_SOURCE_DELIMITER

    def write_header(self):
        """ Write header """
        self.handle.write('#' + '\t'.join(VARIANT_PEPTIDE_TABLE_HEADERS) + '\n')

    def is_valid(self, seq:Seq, canonical_peptides:Set[str],
            cleavage_params:CleavageParams=None):
        """ Check if the peptide is valid """
        min_mw = cleavage_params.min_mw
        min_length = cleavage_params.min_length
        max_length = cleavage_params.max_length
        if SeqUtils.molecular_weight(seq, 'protein') < min_mw:
            return False
        if len(seq) < min_length or len(seq) > max_length:
            return False
        if str(seq) in canonical_peptides:
            return False
        return True

    def add_peptide(self, seq:Seq, peptide_anno:AnnotatedPeptideLabel):
        """ Write """
        start = self.handle.tell()
        for seg in peptide_anno.segments:
            subseq = str(seq[seg.query.start:seg.query.end])
            line = f"{str(seq)}\t{peptide_anno.label}\t{subseq}\t{seg.to_line()}\n"
            self.handle.write(line)
        end = self.handle.tell()
        cur = (start, end)
        if seq in self.index:
            self.index[seq].append(cur)
        else:
            self.index[seq] = [cur]

    def load_peptide(self, seq:Seq):
        """ Load peptide from table """
        labels = set()
        for start, end in self.index[seq]:
            self.handle.seek(start, 0)
            buffer:str = self.handle.read(end - start)
            for line in buffer.rstrip().split('\n'):
                fields = line.split('\t')
                if seq != fields[0]:
                    raise ValueError(
                        f"Peptide ({seq}) do not match with table record ({fields[0]})."
                    )
                labels.add(fields[1])
        label = self.header_delimeter.join(labels)

        return aa.AminoAcidSeqRecord(
            seq=seq,
            description=label,
            name=label
        )

    def load_peptide_annotation(self, seq):
        """ Load peptide annotation """
        anno:Dict[str, AnnotatedPeptideLabel] = {}
        for start, end in self.index[seq]:
            self.handle.seek(start)
            buffer:str = self.handle.read(end - start)
            for line in buffer.rstrip().split('\n'):
                fields = dict(zip(VARIANT_PEPTIDE_TABLE_HEADERS, line.split('\t')))
                if seq != fields['sequence']:
                    raise ValueError(
                        f"Peptide ({seq}) does not match with table record ({fields[0]})."
                    )
                header = fields['header']
                label = anno.setdefault(header, AnnotatedPeptideLabel(label=header, segments=[]))
                query = FeatureLocation(
                    start=int(fields['start']),
                    end=int(fields['end']),
                    start_offset=int(fields['start_offset']),
                    end_offset=int(fields['end_offset'])
                )
                feature_type = fields['feature_type'] if fields['feature_type'] != '.' else None
                feature_id = fields['feature_id'] if fields['feature_id'] != '.' else None

                if fields['ref_start'] != '.':
                    ref_start_raw = int(fields['ref_start'])
                    ref_start = int(ref_start_raw / 3)
                    ref_start_offset = ref_start_raw - ref_start * 3
                    ref_end_raw = int(fields['ref_end'])
                    ref_end = int(ref_end_raw / 3)
                    ref_end_offset = ref_end_raw - ref_end * 3

                    ref = FeatureLocation(
                        start=ref_start,
                        end=ref_end,
                        start_offset=ref_start_offset,
                        end_offset=ref_end_offset
                    )
                else:
                    ref = None
                variant = fields['variant'] if fields['variant'] != '.' else None
                seg = PeptideSegment(
                    query=query, ref=ref, feature_type=feature_type,
                    feature_id=feature_id, variant_id=variant
                )
                label.segments.append(seg)
        return anno

    def write_fasta(self, path:Path):
        """ write fasta """
        with open(path, 'wt') as handle:
            record2title = lambda x: x.description
            writer = FastaIO.FastaWriter(handle, record2title=record2title)
            for seq in self.index:
                peptide = self.load_peptide(seq)
                writer.write_record(peptide)

    def generate_index(self):
        """ Generate index from peptide table file. """
        self.handle.seek(0)
        while True:
            cur_start = self.handle.tell()
            line = self.handle.readline()
            cur_end = cur_start + len(line)
            if not line:
                break
            if line.startswith('#'):
                continue
            seq = line.split('\t')[0]
            indices = self.index.setdefault(seq, [])
            if indices and indices[-1][1] == cur_start:
                indices[-1] = (indices[-1][0], cur_end)
            else:
                indices.append((cur_start, cur_end))
