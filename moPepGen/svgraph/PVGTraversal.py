""" PVGCursor """
from __future__ import annotations
from typing import Tuple, List, Deque, Dict, TYPE_CHECKING, Set
from functools import cmp_to_key
from collections import deque


if TYPE_CHECKING:
    from Bio.Seq import Seq
    from moPepGen.circ.CircRNA import CircRNAModel
    from moPepGen.svgraph.PVGNode import PVGNode
    from moPepGen.svgraph.PVGOrf import PVGOrf
    from moPepGen.seqvar import VariantRecord
    from moPepGen.svgraph.PVGPeptideFinder import PVGPeptideFinder

class PVGCursor():
    """ Helper class for cursors when graph traversal to call peptides. """
    def __init__(self, in_node:PVGNode, out_nodes:Deque[PVGNode], in_cds:bool,
            orfs:List[PVGOrf]=None, cleavage_gain:List[VariantRecord]=None,
            finding_start_site:bool=True):
        """ constructor """
        self.in_node = in_node
        self.out_nodes = out_nodes
        self.in_cds = in_cds
        self.cleavage_gain = cleavage_gain or []
        self.orfs = orfs or []
        self.finding_start_site = finding_start_site

class PVGTraversal():
    """ PVG Traversal. The purpose of this class is to facilitate the graph
    traversal to call variant peptides.
    """
    def __init__(self, check_variants:bool, check_orf:bool,
            pool:PVGPeptideFinder, known_orf_tx:Tuple[int,int]=None,
            known_orf_aa:Tuple[int,int]=None, circ_rna:CircRNAModel=None,
            queue:Deque[PVGCursor]=None,
            stack:Dict[str, Dict[str, PVGCursor]]=None,
            orf_assignment:str='max', backsplicing_only:bool=False,
            find_ass:bool=False,
            reef_kmers:Set[Tuple[str]]=None):
        """ constructor """
        self.check_variants = check_variants
        self.check_orf = check_orf
        self.known_orf_tx = known_orf_tx or (None, None)
        self.known_orf_aa = known_orf_aa or (None, None)
        self.circ_rna = circ_rna
        self.queue = queue or deque([])
        self.pool = pool
        self.stack = stack or {}
        self.orf_assignment = orf_assignment
        self.backsplicing_only = backsplicing_only
        self.find_ass = find_ass
        self.reef_kmers = reef_kmers or set()

    def is_done(self) -> bool:
        """ Check if the traversal is done """
        return not bool(self.queue)

    def has_known_orf(self):
        """ Check if the transcript has any known ORF """
        return self.known_orf_aa[0] is not None

    def known_reading_frame_index(self) -> int:
        """ Get the reading frame index of the known ORF """
        return self.known_orf_tx[0] % 3

    @staticmethod
    def cmp_known_orf_in_frame(x:PVGCursor, y:PVGCursor) -> bool:
        """ comparison for in frame nodes with known ORF """
        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1
        if not x.in_cds and not y.in_cds:
            return -1

        x_start_gain = x.orfs[0].start_gain
        y_start_gain = y.orfs[0].start_gain
        if x_start_gain and not y_start_gain:
            return 1
        if not x_start_gain and y_start_gain:
            return -1
        if x_start_gain and y_start_gain:
            return -1 if sorted(x_start_gain)[0] > sorted(y_start_gain)[0] else 1
        return -1

    @staticmethod
    def cmp_known_orf_frame_shifted(x:PVGCursor, y:PVGCursor) -> bool:
        """ comparison for frameshifing nodes with known ORF """
        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1
        if not x.in_cds and not y.in_cds:
            return -1

        x_start_gain = x.orfs[0].start_gain
        y_start_gain = y.orfs[0].start_gain
        if x_start_gain and not y_start_gain:
            return -1
        if not x_start_gain and y_start_gain:
            return 1
        if x_start_gain and y_start_gain:
            return -1 if sorted(x_start_gain)[0] > sorted(y_start_gain)[0] else 1
        return -1

    @staticmethod
    def cmp_unknown_orf(x:PVGCursor, y:PVGCursor) -> bool:
        """ comparision for unkonwn ORF """
        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1
        if not x.in_cds and not y.in_cds:
            return -1

        x_start_gain = x.orfs[0].start_gain
        y_start_gain = y.orfs[0].start_gain
        if x_start_gain and not y_start_gain:
            return -1
        if not x_start_gain and y_start_gain:
            return 1
        if x_start_gain and y_start_gain:
            return -1 if sorted(x_start_gain)[0] > sorted(y_start_gain)[0] else 1

        if x.cleavage_gain and not y.cleavage_gain:
            return -1
        if not x.cleavage_gain and y.cleavage_gain:
            return 1
        if x.cleavage_gain and y.cleavage_gain:
            return -1 if sorted(x.cleavage_gain)[0] > sorted(y.cleavage_gain)[0] else 1

        return -1

    def cmp_unknown_orf_check_orf(self, x:PVGCursor, y:PVGCursor) -> bool:
        """ comparison for unknown ORF and check ORFs. """
        # pylint: disable=R0911
        if self.find_ass:
            if any(o.orf[0] == self.known_orf_tx[0] for o in x.orfs):
                return -1
            if any(o.orf[0] == self.known_orf_tx[0] for o in y.orfs):
                return 1

        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1
        if not x.in_cds and not y.in_cds:
            return -1

        x_orf = x.orfs[0].orf
        y_orf = y.orfs[0].orf
        if x_orf[0] is not None and (y_orf[0] is None or y_orf[0] == -1):
            return -1
        if (x_orf[0] is None or x_orf[0] == -1) and y_orf[0] is not None:
            return 1
        if x_orf[0] is not None and x_orf[0] != -1 and \
                y_orf[0] is not None and y_orf[0] != -1:
            if x_orf[0] > y_orf[0]:
                return -1 if self.orf_assignment == 'max' else 1
            if x_orf[0] < y_orf[0]:
                return 1 if self.orf_assignment == 'max' else -1

        x_start_gain = x.orfs[0].start_gain
        y_start_gain = y.orfs[0].start_gain
        if x_start_gain and not y_start_gain:
            return -1
        if not x_start_gain and y_start_gain:
            return 1
        if x_start_gain and y_start_gain:
            return -1 if sorted(x_start_gain)[0] > sorted(y_start_gain)[0] else 1
        return -1

    def comp_unknown_orf_keep_all_orfs(self, x:PVGCursor, y:PVGCursor) -> bool:
        """ comparison when all ORFs want to be kept (for circRNA) """
        if self.find_ass:
            if any(o.orf[0] == self.known_orf_tx[0] for o in x.orfs):
                return -1
            if any(o.orf[0] == self.known_orf_tx[0] for o in y.orfs):
                return 1

        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1
        if not x.in_cds and not y.in_cds:
            return -1

        for i, j in zip(x.orfs, y.orfs):
            if i > j:
                return -1
            if i < j:
                return 1

        if len(x.orfs) < len(y.orfs):
            return -1
        if len(x.orfs) > len(y.orfs):
            return 1
        return -1

    def stage(self, in_node:PVGNode, out_nodes:Deque[PVGNode], cursor:PVGCursor):
        """ When a node is visited through a particular edge during the variant
        peptide finding graph traversal, it is staged until all inbond edges
        are visited. """
        in_nodes = self.stack.setdefault(out_nodes[0].id, {})

        if in_node in in_nodes:
            return
        in_nodes[in_node.id] = cursor

        if len(in_nodes) != len(out_nodes[0].in_nodes):
            return

        curs = list(self.stack[out_nodes[0].id].values())

        is_circ_rna = any(x.variant.is_circ_rna() for x in out_nodes[0].variants)

        if self.known_orf_aa[0] is not None:
            if out_nodes[0].reading_frame_index == self.known_reading_frame_index():
                func = self.cmp_known_orf_in_frame
            else:
                func = self.cmp_known_orf_frame_shifted
        elif is_circ_rna:
            func = self.comp_unknown_orf_keep_all_orfs
        elif self.check_orf:
            func = self.cmp_unknown_orf_check_orf
        else:
            func = self.cmp_unknown_orf

        curs.sort(key=cmp_to_key(func))

        cur = curs[0]
        cur.orfs = [x.copy() for x in cur.orfs]
        if is_circ_rna:
            for x in curs[1:]:
                for orf in x.orfs:
                    cur.orfs.append(orf.copy())

        self.queue.appendleft(cur)
