""" Module for peptide variation graph """
from __future__ import annotations
import copy
from typing import Set, Tuple, Deque, Dict
from collections import deque
import itertools
from Bio.Seq import Seq
from moPepGen import aa, svgraph, seqvar, get_equivalent


class PeptideVariantGraph():
    """ Defines the DAG data structure for peptide and variants.

    Attributes:
        root (svgraph.PVGNode): The root node of the graph.
        stop (svgraph.PVGNode): This node is added to the end of each
            branch of the graph to represent the end of the peptide sequence.
            The sequence of this node is '*'
        rule (str): The rule for enzymatic cleavage, e.g., trypsin.
        exception (str): The exception for cleavage rule.
    """
    def __init__(self, root:svgraph.PVGNode):
        """ Construct a PeptideVariantGraph.

        Args:
            root (svgraph.PVGNode): The root node of the graph.
        """
        self.root = root
        self.stop = svgraph.PVGNode(aa.AminoAcidSeqRecord(Seq('*')))
        self.rule = None
        self.exception = None

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

    def cleave_if_posible(self, node:svgraph.PVGNode,
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
            node = node.split_node(site)
            site = node.seq.find_first_enzymatic_cleave_site(
                rule=self.rule,
                exception=self.exception
            )
        return first_node if return_first else node

    def expand_alignment_backward(self, node:svgraph.PVGNode,
            ) -> Tuple[svgraph.PVGNode, Set[svgraph.PVGNode]]:
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
        branches = set()
        to_return = self.stop
        if len(node.in_nodes) > 1:
            raise ValueError('Inbond node has multiple')

        upstream = next(iter(node.in_nodes))

        while node.out_nodes:
            cur = node.out_nodes.pop()
            node.remove_out_edge(cur)
            if cur is self.stop:
                cur = svgraph.PVGNode(
                    seq=node.seq,
                    variants=node.variants,
                    frameshifts=copy.copy(node.frameshifts)
                )
                self.add_stop(cur)
            else:
                cur.seq = node.seq + cur.seq

                if not cur.variants:
                    to_return = cur

                variants = []
                for variant in node.variants:
                    variants.append(variant)
                for variant in cur.variants:
                    variants.append(variant.shift(len(cur.seq)))
                cur.variants = variants

            # only branches out if there is new frameshifts
            frameshifts = copy.copy(node.frameshifts)
            frameshifts_before = len(frameshifts)
            frameshifts.update(cur.frameshifts)
            frameshifts_after = len(frameshifts)

            upstream.add_out_edge(cur)
            last = self.cleave_if_posible(cur)
            if to_return is cur:
                to_return = last

            if frameshifts_before != frameshifts_after:
                branches.add(last)

        upstream.remove_out_edge(node)
        return to_return, branches

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
            cur.remove_out_edge(node)
            cur.add_out_edge(downstream)
            self.cleave_if_posible(node=cur)

        node.remove_out_edge(downstream)
        return downstream

    def merge_join_alignments(self, node:svgraph.PVGNode
            ) -> Tuple[svgraph.PVGNode, Set[svgraph.PVGNode]]:
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
        branches = set()
        to_return = self.stop
        # prev_ref = node.find_reference_prev()
        # if len(prev_ref.in_nodes) > 1:
        #     raise ValueError('upstream has multiple')
        # upstream = next(iter(prev_ref.in_nodes))
        primary_nodes = copy.copy(node.in_nodes)
        secondary_nodes = copy.copy(node.out_nodes)
        for first, second in itertools.product(primary_nodes, secondary_nodes):
            variants = copy.copy(first.variants)
            for variant in second.variants:
                variants.append(variant.shift(len(first.seq) + len(node.seq)))

            if second is self.stop:
                seq = first.seq + node.seq
            else:
                seq = first.seq + node.seq + second.seq

            frameshifts = copy.copy(first.frameshifts)
            frameshifts_before = len(frameshifts)
            frameshifts.update(second.frameshifts)
            frameshifts_after = len(frameshifts)

            new_node = svgraph.PVGNode(
                seq=seq,
                variants=variants,
                frameshifts=frameshifts
            )
            for upstream in first.in_nodes:
                upstream.add_out_edge(new_node)

            if second is self.stop:
                self.add_stop(new_node)
            else:
                for downstream in second.out_nodes:
                    new_node.add_out_edge(downstream)

            if not new_node.variants:
                to_return = new_node

            last_node = self.cleave_if_posible(new_node)

            # return the last cleaved
            if to_return is new_node:
                to_return = last_node

            # aslo add the last cleaved to branch
            if frameshifts_before != frameshifts_after:
                branches.add(new_node)

        for first in primary_nodes:
            self.remove_node(first)

        for second in secondary_nodes:
            self.remove_node(second)

        self.remove_node(node)

        return to_return, branches

    def cross_join_alignments(self, node:svgraph.PVGNode, site:int,
            ) -> Tuple[svgraph.PVGNode, Set[svgraph.PVGNode]]:
        r""" For a given node, split at the given position, expand inbound and
        outbound alignments, and join them with edges.

        In the example below, the node PK is returned

            NCWV           V                 NCWVSTEEK-LPAQV
           /    \         / \               /         X     \
        AER-NCWH-STEEKLPAQ-Q-PK    ->    AER-NCWHSTEEK-LPAQQ-PK

        Args:
            node (svgraph.PVGNode): The node in the graph of target.
        """
        to_return = self.stop
        branches = set()
        left = node
        right = node.split_node(site)
        primary_nodes = copy.copy(left.in_nodes)
        secondary_nodes = copy.copy(right.out_nodes)

        primary_nodes2 = set()
        while primary_nodes:
            first = primary_nodes.pop()
            first.seq = first.seq + left.seq
            while first.out_nodes:
                out_node = first.out_nodes.pop()
                first.remove_out_edge(out_node)
            last = self.cleave_if_posible(first)
            primary_nodes2.add(last)
        primary_nodes = primary_nodes2

        for second in secondary_nodes:
            second.seq = right.seq + second.seq

            variants = []
            for variant in right.variants:
                variants.append(variant)
            for variant in second.variants:
                variants.append(variant.shift(len(right.seq)))
            second.variants = variants

            while second.in_nodes:
                in_node = second.in_nodes.pop()
                in_node.remove_out_edge(second)

            # only branch out when there is new frameshifts
            for frameshift in second.frameshifts:
                if frameshift not in right.frameshifts:
                    branches.add(second)

            last = self.cleave_if_posible(second)
            if not variants:
                to_return = last

        self.remove_node(left)
        self.remove_node(right)

        for first, second in itertools.product(primary_nodes, secondary_nodes):
            first.add_out_edge(second)

        return to_return, branches

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

            # each is a branch
            if cur.seq is None:
                for out_node in cur.out_nodes:
                    queue.appendleft(out_node)
                continue

            site = cur.seq.find_first_enzymatic_cleave_site(rule=self.rule,
                exception=self.exception)
            while site > -1:
                cur = cur.split_node(site)
                site = cur.seq.find_first_enzymatic_cleave_site(rule=self.rule,
                    exception=self.exception)

            if self.next_is_stop(cur):
                if len(cur.out_nodes) > 1:
                    # return branch
                    cur, branches = self.expand_alignment_backward(cur)
                    for branch in branches:
                        queue.appendleft(branch)
                continue

            # expand the variant alignments to the left
            # returns the new reference node in the alignment.
            # This is when the node between to bubbles don't have any cleavages
            if len(cur.in_nodes) == 1:
                cur, branches = self.expand_alignment_backward(cur)
                for branch in branches:
                    queue.appendleft(branch)

                if cur is self.stop:
                    continue

                if self.next_is_stop(cur):
                    # leave it to the next iteration to handle
                    queue.append(cur)
                    continue
                cur = cur.find_reference_next()

            if cur is None:
                raise ValueError()
            sites = cur.seq.find_all_enzymatic_cleave_sites(rule=self.rule,
                exception=self.exception)

            branches = set()
            if len(sites) == 0:
                if self.next_is_stop(cur):
                    self.expand_alignment_forward(cur)
                    continue
                cur, branches = self.merge_join_alignments(cur)
                if len(cur.out_nodes) == 1 and not self.next_is_stop(cur):
                    cur = cur.find_reference_next()
            elif len(sites) == 1:
                if self.next_is_stop(cur):
                    cur.split_node(sites[0])
                    cur = self.expand_alignment_forward(cur)
                    queue.appendleft(cur)
                    continue
                cur, branches = self.cross_join_alignments(cur, sites[0])
                if not self.next_is_stop(cur):
                    cur = cur.find_reference_next()
            else:
                cur.split_node(sites[0])
                cur = self.expand_alignment_forward(cur)
            for branch in branches:
                queue.appendleft(branch)
            if cur is not self.stop:
                queue.appendleft(cur)

    def call_vaiant_peptides(self, miscleavage:int=2
            ) -> Set[aa.AminoAcidSeqRecord]:
        """ Walk through the graph and find all variated peptides.

        Args:
            miscleavage (int): Number of miscleavages allowed.

        Return:
            A set of aa.AminoAcidSeqRecord.
        """
        queue:Deque[svgraph.PVGNode] = deque([self.root])
        visited:Set[svgraph.PVGNode] = set()
        variant_peptides:Set[aa.AminoAcidSeqRecord] = set()
        peptide_ids:Dict[str, int] = {}
        while queue:
            cur = queue.pop()
            if cur.seq is None:
                for out_node in cur.out_nodes:
                    queue.appendleft(out_node)
                continue

            # skip the node if already visited
            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                continue

            # add out nodes to the queue
            for out_node in cur.out_nodes:
                if out_node is not self.stop:
                    queue.appendleft(out_node)

            node_paths:Deque[Deque[svgraph.PVGNode]]
            node_paths = deque([deque([cur])])
            i = 0
            while i < miscleavage and node_paths:
                i += 1
                next_batch:Deque[Deque[svgraph.PVGNode]] = deque()
                while node_paths:
                    nodes = node_paths.pop()

                    variants = []
                    seq = None
                    for node in nodes:
                        if seq is None:
                            seq = node.seq
                        else:
                            seq = seq + node.seq
                        if node.variants:
                            variants.extend(node.variants)

                    if not variants:
                        continue

                    variant_label = ''
                    variant:seqvar.VariantRecordWithCoordinate
                    for variant in variants:
                        variant_label += ('|' + str(variant.variant.id))

                    seq.id = seq.transcript_id
                    seq.name = seq.transcript_id
                    seq.description = seq.transcript_id + variant_label

                    if seq.description not in peptide_ids:
                        peptide_ids[seq.description] = 0

                    peptide_ids[seq.description] += 1
                    seq.description += '|' + str(peptide_ids[seq.description])

                    same_peptide = get_equivalent(variant_peptides, seq)
                    if same_peptide:
                        same_peptide:Seq
                        same_peptide.description += '||' + seq.description
                    else:
                        variant_peptides.add(seq)

                    if i + 1 < miscleavage:
                    # create a new batch of node paths for the next iteration
                        last_node = nodes[-1]
                        for node in last_node.out_nodes:
                            if node is self.stop:
                                continue
                            new_list = copy.copy(nodes)
                            new_list.append(node)
                            next_batch.appendleft(new_list)

                node_paths = next_batch

        return variant_peptides
