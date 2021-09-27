""" Module for peptide variation graph """
from __future__ import annotations
import copy
from typing import Set, Deque, Dict, List, Tuple
from collections import deque
import itertools
from Bio.Seq import Seq
from moPepGen import aa, svgraph
from moPepGen.svgraph.VariantPeptideDict import VariantPeptideDict


class PeptideVariantGraph():
    """ Defines the DAG data structure for peptide and variants.

    Attributes:
        root (svgraph.PVGNode): The root node of the graph.
        stop (svgraph.PVGNode): This node is added to the end of each
            branch of the graph to represent the end of the peptide sequence.
            The sequence of this node is '*'
        id (str): The ID of the transcript. For circRNA, this could be the
            circRNA ID.
        has_known_orf (bool): Whether the given transcript has any known
            ORF. If this is a noncoding transcript, 3-frame searching is
            performed.
        rule (str): The rule for enzymatic cleavage, e.g., trypsin.
        exception (str): The exception for cleavage rule.
        orfs (set): A set of all unique ORF start positions.
        reading_frames (List[PVGNode]): The 3 reading frames, each being the
            start node or None if not avaialable.
        orf_id_map (Dict[int, str]): An ID map of ORF IDs from ORF start
            position.
    """
    def __init__(self, root:svgraph.PVGNode, _id:str,
            known_orf:List[int,int], rule:str=None, exception:str=None,
            orfs:Set[int]=None, reading_frames:List[svgraph.PVGNode]=None,
            orf_id_map:Dict[int,str]=None):
        """ Construct a PeptideVariantGraph """
        self.root = root
        self.id = _id
        self.known_orf = known_orf
        self.stop = svgraph.PVGNode(aa.AminoAcidSeqRecord(Seq('*')), None)
        self.rule = rule
        self.exception = exception
        self.orfs = orfs or set()
        self.reading_frames = reading_frames or [None, None, None]
        self.orf_id_map = orf_id_map or {}

    def add_stop(self, node:svgraph.PVGNode):
        """ Add the stop node after the specified node. """
        node.add_out_edge(self.stop)

    def next_is_stop(self, node:svgraph.PVGNode) -> bool:
        """ Checks if stop is linked as a outbound node. """
        return self.stop in node.out_nodes

    @staticmethod
    def remove_node(node:svgraph.PVGNode) -> None:
        """ Remove a given node, and also remove the all edges from and to it.
        """
        while node.in_nodes:
            in_node = node.in_nodes.pop()
            in_node.remove_out_edge(node)
        while node.out_nodes:
            out_node = node.out_nodes.pop()
            node.remove_out_edge(out_node)
        del node

    def known_reading_frame_index(self) -> int:
        """ """
        return self.known_orf[0] % 3

    def has_known_orf(self) -> bool:
        """ Check if it has a known ORF """
        return self.known_orf[0]

    def cleave_if_possible(self, node:svgraph.PVGNode,
            return_first:bool=False) -> svgraph.PVGNode:
        """ For a given node, it checks whether there is any cleave site and
        split at each site.

        Args:
            node (svgraph.PVGNode): The target node to cleave
            return_first (bool): When true, returns the first node rather
                than last. Defaults to False
        """
        first_node = node
        site = node.seq.find_first_enzymatic_cleave_site(
            rule=self.rule,
            exception=self.exception
        )
        while site > -1:
            node = node.split_node(site, cleavage=True)
            site = node.seq.find_first_enzymatic_cleave_site(
                rule=self.rule,
                exception=self.exception
            )
        return first_node if return_first else node

    def group_outbond_siblings(self, upstream:svgraph.PVGNode
            ) -> List[List[svgraph.PVGNode]]:
        """ Group the downstream nodes with variant alignments if they are
        connected to the same downstream node. Sibling nodes are the ones that
        share the same upstream and downstream nodes.

        Args:
            upstream (svgraph.TVGNode): The upstream node of which the outbond
                nodes are grouped.

        Return:
            A 2-dimensional list, that each child list contains nodes that are
            sibling to each other (ie sharing the same upstream and downstream
            node).
        """
        out_nodes = copy.copy(upstream.out_nodes)
        groups:List[List[svgraph.PVGNode]] = []
        while out_nodes:
            out_node = out_nodes.pop()
            if out_node is self.stop:
                continue
            downstream = out_node.find_reference_next()
            group:List[svgraph.TVGNode] = [out_node]
            if downstream and downstream is not self.stop:
                for node in downstream.in_nodes:
                    if node in out_nodes:
                        group.append(node)
                        out_nodes.remove(node)
            groups.append(group)
        return groups

    def expand_alignment_backward(self, node:svgraph.PVGNode,
            ) -> Set[svgraph.PVGNode]:
        r""" Expand the variant alignment bubble backward to the previous
        cleave site. The sequence of the input node is first prepended to each
        of outbound node, and then the inbound node of those outbond nodes are
        then pointed to the inbond node of the input node.

        In the example below, the node H and V are expanded backward to include
        NCW. The input variable node should be the node NCW. The returned node
        is ST

                V                    NCWV
               / \                  /    \
        AER-NCW-H-ST     ->      AER-NCWH-ST

        Args:
            node (svgraph.Peptide): The node that the outbond variant alignment
            bubble should be expanded.

        Returns:
            A tuple of size 2 is returned. The first element is the end node
            of the variant alignment bubble, and the second is a set of all
            branches with frameshifting variants.
        """
        upstream = next(iter(node.in_nodes))

        downstreams = set()

        for siblings in self.group_outbond_siblings(node):
            for cur in siblings:
                node.remove_out_edge(cur)
                if cur is self.stop:
                    cur = svgraph.PVGNode(
                        seq=copy.copy(node.seq),
                        variants=copy.copy(node.variants),
                        frameshifts=copy.copy(node.frameshifts),
                        orf=node.orf,
                        reading_frame_index=node.reading_frame_index
                    )
                    self.add_stop(cur)
                else:
                    cur.seq = node.seq + cur.seq

                    variants = copy.copy(node.variants)
                    for variant in cur.variants:
                        variants.append(variant.shift(len(node.seq)))
                    cur.variants = variants

                frameshifts = copy.copy(node.frameshifts)
                frameshifts.update(cur.variants)
                cur.frameshifts = frameshifts

                upstream.add_out_edge(cur)
                last = self.cleave_if_possible(cur)

            if len(siblings) == 1:
                downstreams.add(last)
            else:
                if len(last.out_nodes) > 1:
                    raise ValueError('Multiple out nodes found.')
                downstreams.add(list(last.out_nodes)[0])

        upstream.remove_out_edge(node)
        # return to_return, branches
        return downstreams

    def expand_alignment_forward(self, node:svgraph.PVGNode
            ) -> svgraph.PVGNode:
        r""" Expand the upsteam variant alignment bubble to the end of the
        node of input. The sequencing of the node of input is first appended
        to each of the leading nodes, and the outbound node of those nodes
        are then pointed to the input node's outbound node.

        For each of the expanded node, cleavage sites are searched. They will
        be spliced if any.

        In the example below, nodes NCWV and NCWH are expended to include
        STQQPK. The node QQPQE is then returned.

            NCWV                           NCWVSTQQPK
           /    \                         /          \
        AER-NCWH-STQQPK-QQPQE    ->    AER-NCWHSTQQPK-QQPQE

        Args:
            node (svgraph.PVGNode): The node that the inbound variant
                alignment bubble to be expanded.

        Returns:
            The end node is returned.
        """
        if len(node.out_nodes) > 1:
            raise ValueError('Outbond node has multiple')

        downstream = next(iter(node.out_nodes))

        while node.in_nodes:
            cur = node.in_nodes.pop()
            cur.seq = cur.seq + node.seq

            for variant in node.variants:
                cur.variants.append(variant.shift(len(cur.seq)))
            cur.frameshifts.update(node.variants)

            cur.remove_out_edge(node)
            cur.add_out_edge(downstream)
            self.cleave_if_possible(node=cur)

        node.remove_out_edge(downstream)
        return downstream

    def merge_join_alignments(self, node:svgraph.PVGNode
            ) -> Set[svgraph.PVGNode]:
        r""" For a given node, join all the inbond and outbond nodes with any
        combinations.

        In the example below, node NCWHSTQQ is returned

                                       NCWHSTQV
                                      /        \
            NCWV     V               | NCWVSTQV |
           /    \   / \              |/        \|
        AER-NCWH-STQ-Q-PK    ->    AER-NCWHSTQQ-PK
                                      \        /
                                       NCWVSTQQ

        Args:
            node (svgraph.PVGNode): The node in the graph of target.
         """
        downstreams = set()

        primary_nodes = copy.copy(node.in_nodes)
        trash = {node}

        for siblings in self.group_outbond_siblings(node):
            if len(siblings) == 1:
                sibling = siblings[0]
                node.remove_out_edge(sibling)
                for first in primary_nodes:
                    seq = first.seq + node.seq
                    variants = copy.copy(first.variants)
                    for variant in node.variants:
                        variants.append(variant.shift(len(first.seq)))
                    frameshifts = copy.copy(first.frameshifts)
                    frameshifts.update(node.variants)

                    new_node = svgraph.PVGNode(
                        seq=seq, variants=variants, frameshifts=frameshifts,
                        orf=first.orf, reading_frame_index=first.reading_frame_index
                    )
                    for upstream in first.in_nodes:
                        upstream.add_out_edge(new_node)
                    new_node.add_out_edge(sibling)
                    trash.add(first)
                downstreams.add(sibling)
            else:
                for first, second in itertools.product(primary_nodes, siblings):

                    variants = copy.copy(first.variants)
                    for variant in node.variants:
                        variants.append(variant.shift(len(first.seq)))
                    for variant in second.variants:
                        variants.append(variant.shift(len(first.seq) + len(node.seq)))

                    if second is self.stop:
                        seq = first.seq + node.seq
                    else:
                        seq = first.seq + node.seq + second.seq

                    frameshifts = copy.copy(first.frameshifts)
                    frameshifts.update(node.variants)
                    frameshifts.update(second.variants)

                    new_node = svgraph.PVGNode(
                        seq=seq,
                        variants=variants,
                        frameshifts=frameshifts,
                        orf=first.orf,
                        reading_frame_index=first.reading_frame_index
                    )
                    for upstream in first.in_nodes:
                        upstream.add_out_edge(new_node)

                    if second is self.stop:
                        self.add_stop(new_node)
                    else:
                        for downstream in second.out_nodes:
                            new_node.add_out_edge(downstream)

                    last = self.cleave_if_possible(new_node)

                    trash.add(first)
                    trash.add(second)

                if len(last.out_nodes) > 1:
                    raise ValueError('Multiple out nodes found.')
                downstreams.add(list(last.out_nodes)[0])

        for x in trash:
            self.remove_node(x)

        return downstreams

    def cross_join_alignments(self, node:svgraph.PVGNode, site:int,
            ) -> List[svgraph.PVGNode]:
        r""" For a given node, split at the given position, expand inbound and
        outbound alignments, and join them with edges.

        In the example below, the node PK is returned

            NCWV           V                 NCWVSTEEK-LPAQV
           /    \         / \               /         X     \
        AER-NCWH-STEEKLPAQ-Q-PK    ->    AER-NCWHSTEEK-LPAQQ-PK

        Args:
            node (svgraph.PVGNode): The node in the graph of target.
        """
        left = node
        right = node.split_node(site, cleavage=True)
        primary_nodes = copy.copy(left.in_nodes)
        secondary_nodes = copy.copy(right.out_nodes)

        downstreams = set()

        primary_nodes2 = set()
        while primary_nodes:
            first = primary_nodes.pop()
            first.seq = first.seq + left.seq

            for variant in left.variants:
                first.variants.append(variant.shift(len(first.seq)))
            first.frameshifts.update(left.variants)

            while first.out_nodes:
                out_node = first.out_nodes.pop()
                first.remove_out_edge(out_node)
            last = self.cleave_if_possible(first)
            primary_nodes2.add(last)
        primary_nodes = primary_nodes2

        for siblings in self.group_outbond_siblings(right):
            for second in siblings:
                second.seq = right.seq + second.seq

                variants = copy.copy(right.variants)
                for variant in second.variants:
                    variants.append(variant.shift(len(right.seq)))
                second.variants = variants
                second.frameshifts.update(right.variants)

                while second.in_nodes:
                    in_node = second.in_nodes.pop()
                    in_node.remove_out_edge(second)

                second.cleavage = True

                last = self.cleave_if_possible(second)

            if len(siblings) == 1:
                downstreams.add(last)
            else:
                if len(last.out_nodes) > 1:
                    raise ValueError('Multiple out nodes found.')
                downstreams.add(list(last.out_nodes)[0])

        self.remove_node(left)
        self.remove_node(right)

        for first, second in itertools.product(primary_nodes, secondary_nodes):
            first.add_out_edge(second)

        return downstreams

    def form_cleavage_graph(self, rule:str, exception:str=None) -> None:
        """ Form a cleavage graph from a variant graph. After calling this
        method, every each in the graph should represent a cleavage in the
        sequence of the pull peptide of reference and variated.

        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
        """
        self.rule = rule
        self.exception = exception
        queue = deque([self.root])
        branch = set()

        while queue:
            cur = queue.pop()

            if cur is self.stop:
                continue

            if cur.seq is None:
                for out_node in cur.out_nodes:
                    queue.appendleft(out_node)
                continue

            if len(cur.in_nodes) == 1:

                site = cur.seq.find_first_enzymatic_cleave_site(rule=self.rule,
                    exception=self.exception)
                while site > -1:
                    cur = cur.split_node(site, cleavage=True)
                    site = cur.seq.find_first_enzymatic_cleave_site(
                        rule=self.rule, exception=self.exception)

                if self.next_is_stop(cur):
                    if len(cur.out_nodes) > 1:
                        # return branch
                        branches = self.expand_alignment_backward(cur)
                        for branch in branches:
                            queue.appendleft(branch)
                    continue

                branches = self.expand_alignment_backward(cur)
                for branch in branches:
                    if not branch.is_bridge():
                        queue.appendleft(branch)
                continue

            sites = cur.seq.find_all_enzymatic_cleave_sites(rule=self.rule,
                exception=self.exception)

            if len(sites) == 0:
                if self.next_is_stop(cur):
                    self.expand_alignment_forward(cur)
                    continue
                if cur.cleavage:
                    branches = self.expand_alignment_backward(cur)
                else:
                    branches = self.merge_join_alignments(cur)
                for branch in branches:
                    if not branch.is_bridge():
                        queue.appendleft(branch)
            elif len(sites) == 1:
                if self.next_is_stop(cur):
                    cur.split_node(sites[0], cleavage=True)
                    cur = self.expand_alignment_forward(cur)
                    queue.appendleft(cur)
                    continue
                branches = self.cross_join_alignments(cur, sites[0])
                for branch in branches:
                    if not branch.is_bridge():
                        queue.appendleft(branch)
            else:
                right = cur.split_node(sites[0], cleavage=True)
                cur = self.expand_alignment_forward(cur)
                site = right.seq.find_first_enzymatic_cleave_site(
                    rule=self.rule, exception=self.exception)
                right = right.split_node(site, cleavage=True)
                queue.appendleft(right)

    def init_orfs(self) -> None:
        """ Init the orfs list with all existing orfs. """
        self.orfs = {node.orf[0] for node in self.root.out_nodes}

    def update_orf(self, node:svgraph.PVGNode) -> None:
        """ Update the orf list with the orf start position of the given node.
        """
        self.orfs.add(node.orf[0])

    def create_orf_id_map(self) -> None:
        """ Creates and map for ORF start site and its ID (ORF1, ORF2, etc) """
        orfs = list(self.orfs)
        orfs.sort()
        self.orf_id_map = {v:f"ORF{i+1}" for i,v in enumerate(orfs)}


    def call_variant_peptides(self, miscleavage:int=2, check_variants:bool=True,
            check_orf:bool=False, keep_all_occurrence:bool=True
            ) -> Set[aa.AminoAcidSeqRecord]:
        """ Walk through the graph and find all variated peptides.

        Args:
            miscleavage (int): Number of miscleavages allowed.
            check_variants (bool): When true, only peptides that carries at
                least 1 variant are kept. And when false, all unique peptides
                are reported.
            check_orf (bool): When true, the ORF ID will be added to the
                variant peptide label.
            keep_all_occurrence: Whether to keep all occurance of the peptide
                within the graph/transcript.

        Return:
            A set of aa.AminoAcidSeqRecord.
        """
        self.init_orfs()

        queue:Deque[Tuple[svgraph.PVGNode,bool]] = deque([(self.root,True)])
        visited:Set[svgraph.PVGNode] = set()
        peptide_pool = VariantPeptideDict(tx_id=self.id)

        if self.has_known_orf():
            known_orf_aa = (
                int((self.known_orf[0] - self.known_reading_frame_index())/3),
                int((self.known_orf[1] - self.known_reading_frame_index())/3)
            )

        while queue:
            cur, in_cds = queue.pop()
            orf = cur.orf
            if cur is self.root and cur.seq is None:
                if self.has_known_orf():
                    i = self.known_reading_frame_index()
                    out_node = self.reading_frames[i]
                    node_start = out_node.seq.locations[0].ref.start
                    in_cds = node_start >= known_orf_aa[0]
                    queue.appendleft((out_node, in_cds))
                else:
                    for out_node in cur.out_nodes:
                        queue.appendleft((out_node,False))
                continue

            if cur is self.stop:
                continue

            # skip the node if already visited
            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                continue

            if self.has_known_orf():
                # add out nodes to the queue
                if in_cds:
                    stop_index = cur.seq.get_query_index(known_orf_aa[1])
                    if stop_index > -1:
                        cur.truncate_right(stop_index)

                    truncate_index = cur.seq.seq.find('*')
                    if truncate_index > -1:
                        cur.truncate_right(truncate_index)
                        cur.truncated = True

                    peptide_pool.add_miscleaved_sequences(cur, miscleavage,
                        check_variants)
                    if stop_index == -1 and truncate_index == -1:
                        for out_node in cur.out_nodes:
                            if out_node is not self.stop:
                                queue.appendleft((out_node, True))
                else:
                    start_index = cur.seq.get_query_index(known_orf_aa[0])
                    if start_index == -1:
                        out_node = cur.find_reference_next()
                        queue.appendleft((out_node,False))
                    else:
                        visited.remove(cur)
                        cur.truncate_left(start_index)
                        queue.appendleft((cur, True))
                continue

            node_list = []
            trash = set()

            orf = cur.orf
            stop_index = cur.seq.seq.find('*')
            start_index = cur.seq.seq.find('M')
            if in_cds and stop_index == -1:
                node_list.append(cur)
            if in_cds and stop_index > -1:
                cur_copy = cur.copy()
                cur_copy.truncate_right(stop_index)
                cur_copy.remove_out_edges()
                node_list.append(cur_copy)
                trash.add(cur_copy)
                if stop_index > start_index:
                    in_cds = False
                    orf = [None, None]
            last_stop_index, last_start_index = None, None
            while stop_index > -1 or start_index > -1:
                if 0 < start_index < stop_index and \
                        (start_index != last_start_index or \
                        stop_index != last_stop_index):
                    cur_copy = cur[start_index:stop_index]
                    cur_copy.remove_out_edges()
                    orf_start = cur_copy.get_orf_start()
                    cur_copy.orf = [orf_start, None]
                    self.add_stop(cur_copy)
                    self.update_orf(cur_copy)
                    node_list.append(cur_copy)
                    trash.add(cur_copy)
                if start_index > -1 and stop_index == -1 \
                        and start_index != last_start_index:
                    cur_copy = cur.copy()
                    cur_copy.truncate_left(start_index)
                    orf_start = cur_copy.get_orf_start()
                    cur_copy.orf = [orf_start, None]
                    self.update_orf(cur_copy)
                    if not in_cds:
                        in_cds = True
                        orf = cur_copy.orf
                    node_list.append(cur_copy)
                    trash.add(cur_copy)

                last_start_index, last_stop_index = start_index, stop_index

                if -1 < start_index < stop_index or (stop_index == -1 \
                        and start_index > -1):
                    x = cur.seq.seq[start_index+1:].find('M')
                    start_index = x + start_index + 1 if x > -1 else x
                elif -1 < stop_index < start_index or (stop_index > -1 \
                        and start_index == -1):
                    x = cur.seq.seq[stop_index+1:].find('*')
                    stop_index = x + stop_index + 1 if x > -1 else x

            for out_node in cur.out_nodes:
                if out_node is not self.stop:
                    out_node.orf = orf
                    queue.appendleft((out_node, in_cds))

            for node in node_list:
                peptide_pool.add_miscleaved_sequences(node, miscleavage,
                    check_variants)
            for node in trash:
                self.remove_node(node)

        if check_orf:
            self.create_orf_id_map()

        return peptide_pool.get_peptide_sequences(
            keep_all_occurrence=keep_all_occurrence, orf_id_map=self.orf_id_map
        )
