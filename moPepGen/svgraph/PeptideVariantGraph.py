""" Module for peptide variation graph """
from __future__ import annotations
import copy
from typing import Iterable, Set, Deque, Dict, List, Tuple
from collections import deque
from functools import cmp_to_key
import itertools
from Bio.Seq import Seq
from moPepGen import aa, seqvar
from moPepGen.svgraph.VariantPeptideDict import VariantPeptideDict
from moPepGen.svgraph.PVGNode import PVGNode


T = Tuple[Set[PVGNode],Dict[PVGNode,List[PVGNode]]]

class PeptideVariantGraph():
    """ Defines the DAG data structure for peptide and variants.

    Attributes:
        root (PVGNode): The root node of the graph.
        stop (PVGNode): This node is added to the end of each
            branch of the graph to represent the end of the peptide sequence.
            The sequence of this node is '*'
        id (str): The ID of the transcript. For circRNA, this could be the
            circRNA ID.
        has_known_orf (bool): Whether the givenf transcript has any known
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
    def __init__(self, root:PVGNode, _id:str,
            known_orf:List[int,int], rule:str=None, exception:str=None,
            orfs:Set[Tuple[int,int]]=None, reading_frames:List[PVGNode]=None,
            orf_id_map:Dict[int,str]=None):
        """ Construct a PeptideVariantGraph """
        self.root = root
        self.id = _id
        self.known_orf = known_orf
        self.stop = PVGNode(aa.AminoAcidSeqRecord(Seq('*')), None)
        self.rule = rule
        self.exception = exception
        self.orfs = orfs or set()
        self.reading_frames = reading_frames or [None, None, None]
        self.orf_id_map = orf_id_map or {}

    def add_stop(self, node:PVGNode):
        """ Add the stop node after the specified node. """
        node.add_out_edge(self.stop)

    def next_is_stop(self, node:PVGNode) -> bool:
        """ Checks if stop is linked as a outbound node. """
        return self.stop in node.out_nodes

    @staticmethod
    def remove_node(node:PVGNode) -> None:
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
        """ Get the readign frame index of the known ORF. """
        return self.known_orf[0] % 3

    def has_known_orf(self) -> bool:
        """ Check if it has a known ORF """
        return self.known_orf[0] is not None

    def cleave_if_possible(self, node:PVGNode,
            return_first:bool=False) -> PVGNode:
        """ For a given node, it checks whether there is any cleave site and
        split at each site.

        Args:
            node (PVGNode): The target node to cleave
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

    @staticmethod
    def find_routes_for_merging(node:PVGNode, cleavage:bool=False
            ) -> Set[Tuple[PVGNode]]:
        """ Find all start and end nodes for merging.

        Args:
            node (PVGNode)
            cleavage (bool): When True, the input node is treated as the start
                node, so it's in-bond nodes are not looked.
        """
        routes:Set[Tuple[PVGNode]] = set()
        visited:Set[PVGNode] = {node}
        if cleavage:
            for out_node in node.out_nodes:
                route = (node, out_node)
                routes.add(route)
                visited.add(out_node)
        else:
            for in_node, out_node in itertools.product(node.in_nodes, node.out_nodes):
                route = (in_node, node, out_node)
                routes.add(route)
                visited.add(in_node)
                visited.add(out_node)

        for cur in copy.copy(routes):
            for in_node in cur[-1].in_nodes:
                if in_node.is_bridge() and in_node not in visited:
                    for out_node in cur[-1].out_nodes:
                        new_route = (in_node, cur[-1],)
                        routes.add(new_route)
                        visited.add(in_node)
        return routes

    def merge_nodes_routes(self, routes:Set[Tuple[PVGNode]]
            ) -> Tuple[Set[PVGNode], Dict[PVGNode, List[PVGNode]]]:
        """ For any nodes between given serious of start and end nodes, merge
        them if they are connected. The total number of resulted merged nodes
        equals to the number of nodes in the end set. This function is used
        primary in creating the cleavage graph.

        If a merged node is a singleton (i.e., no other node sharing the same
        end node with it in the result), itself is returned, otherwise, returns
        the next node.

        Args:
            start_nodes (Iterable[PVGNode]): serious of nodes to start merging.
            end_nodes (Iterable[PVGNode]): serious of nodes to merge to.
        Return:
            A set of nodes to be used in the next iteration of creating the
            cleavage graph.
        """
        trash:Set[Tuple[PVGNode,PVGNode]] = set()
        new_nodes:Set[PVGNode] = set()
        start_nodes = {route[0] for route in routes}
        end_nodes = {route[-1] for route in routes}
        inbridge_list:Dict[PVGNode,List[PVGNode]] = {}

        # for end_node in end_nodes:
        #     site = end_node.seq.find_first_enzymatic_cleave_site(self.rule)
        #     if site > -1:
        #         end_node.split_node(site, True)

        for route in routes:
            for i,node in enumerate(route):
                if i == 0:
                    new_node = node.copy(in_nodes=False, out_nodes=False)
                else:
                    new_node.append_right(node)
                    trash.add((route[i-1], node))

            for in_node in route[0].in_nodes:
                in_node.add_out_edge(new_node)

            for out_node in route[-1].out_nodes:
                new_node.add_out_edge(out_node)

            new_nodes.add(new_node)

            if route[0].is_bridge():
                val = inbridge_list.setdefault(route[0], [])
                val.append(new_node)

        for in_node,out_node in trash:
            in_node.remove_out_edge(out_node)

        for node in start_nodes:
            if not node.out_nodes:
                self.remove_node(node)

        for node in end_nodes:
            if not node.in_nodes:
                self.remove_node(node)

        for new_node in copy.copy(new_nodes):
            last = self.cleave_if_possible(new_node)
            new_nodes.remove(new_node)
            new_nodes.add(last)

        return new_nodes, inbridge_list

    def move_downstreams(self, nodes:Iterable[PVGNode], reading_frame_index:int
            ) -> Set[PVGNode]:
        """ Get the downstream nodes to be added to the queue for creating
        cleavage graph. This function should be called after the
        `merge_nodes_routes` in the graph updating functions such as
        merge_forward, merge_backward etc. """
        downstreams:Set[PVGNode] = set()
        for node in nodes:
            if node.reading_frame_index != reading_frame_index:
                continue
            if node.is_bridge():
                continue
            is_sharing_downstream = any(any(y is not node and y in nodes \
                for y in x.in_nodes) for x in node.out_nodes if x is not self.stop)
            if not is_sharing_downstream:
                downstreams.add(node)
            else:
                for downstream in node.out_nodes:
                    if downstream is self.stop:
                        continue
                    downstreams.add(downstream)
        return downstreams

    def expand_backward(self, node:PVGNode) -> T:
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
            node (PVGNode): The node that the outbond variant alignment
            bubble should be expanded.

        Returns:
            A set of nodes that are either the resulted merged node if it is
            singleton, or the downstream node otherwise.
        """
        reading_frame_index = node.reading_frame_index
        routes = self.find_routes_for_merging(node, True)
        new_nodes, inbridge_list = self.merge_nodes_routes(routes)
        downstreams = self.move_downstreams(new_nodes, reading_frame_index)
        return downstreams, inbridge_list

    def expand_forward(self, node:PVGNode) -> T:
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
            node (PVGNode): The node that the inbound variant
                alignment bubble to be expanded.

        Returns:
            The end node is returned.
        """
        routes = set()
        reading_frame_index = node.reading_frame_index
        for in_node in node.in_nodes:
            routes.add((in_node, node))
        new_nodes, inbridge_list = self.merge_nodes_routes(routes)
        downstreams = self.move_downstreams(new_nodes, reading_frame_index)
        return downstreams, inbridge_list

    def merge_join(self, node:PVGNode) -> T:
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
            node (PVGNode): The node in the graph of target.
        returns:
            A set of nodes that are either the resulted merged node if it is
            singleton, or the downstream node otherwise.
        """
        reading_frame_index = node.reading_frame_index
        routes = self.find_routes_for_merging(node)
        new_nodes, inbridge_list = self.merge_nodes_routes(routes)
        downstreams = self.move_downstreams(new_nodes, reading_frame_index)
        return downstreams, inbridge_list

    def cross_join(self, node:PVGNode, site:int) -> T:
        r""" For a given node, split at the given position, expand inbound and
        outbound alignments, and join them with edges.

        In the example below, the node PK is returned

            NCWV           V                 NCWVSTEEK-LPAQV
           /    \         / \               /         X     \
        AER-NCWH-STEEKLPAQ-Q-PK    ->    AER-NCWHSTEEK-LPAQQ-PK

        Args:
            node (PVGNode): The node in the graph of target.
        """
        reading_frame_index = node.reading_frame_index
        head = node.split_node(site, cleavage=True)

        # there shouldn't have any bridge out in inbond nodes
        self.expand_forward(node)

        routes = self.find_routes_for_merging(head, True)
        new_nodes, inbridge_list = self.merge_nodes_routes(routes)
        downstreams = self.move_downstreams(new_nodes, reading_frame_index)
        return downstreams, inbridge_list

    def create_cleavage_graph(self, rule:str, exception:str=None) -> None:
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

        inbridge_list:Dict[PVGNode,List[PVGNode]] = {}

        while queue:
            cur = queue.pop()

            if cur is self.stop:
                continue

            if cur.seq is None:
                for out_node in cur.out_nodes:
                    queue.appendleft(out_node)
                continue

            if cur in inbridge_list:
                for node in inbridge_list[cur]:
                    queue.appendleft(node)
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
                        branches, inbridges = self.expand_backward(cur)
                        for branch in branches:
                            queue.appendleft(branch)
                        for key, val in inbridges.items():
                            inbridge_list[key] = val
                    continue

                if cur.is_bridge():
                    continue
                if len(cur.out_nodes) == 1:
                    downstream = list(cur.out_nodes)[0]
                    site = downstream.seq.find_first_enzymatic_cleave_site(
                        self.rule, self.exception
                    )
                    if site > -1:
                        downstream.split_node(site, True)
                    branches, inbridges = self.expand_forward(downstream)
                else:
                    branches, inbridges = self.expand_backward(cur)
                for branch in branches:
                    if cur.reading_frame_index == branch.reading_frame_index and \
                            not branch.is_bridge():
                        queue.appendleft(branch)
                for key, val in inbridges.items():
                    inbridge_list[key] = val
                continue

            sites = cur.seq.find_all_enzymatic_cleave_sites(rule=self.rule,
                exception=self.exception)

            if len(sites) == 0:
                if self.next_is_stop(cur) or cur.is_bridge():
                    if not cur.cleavage:
                        _,inbridges = self.expand_forward(cur)
                        for key, val in inbridges.items():
                            inbridge_list[key] = val
                    continue
                if cur.cleavage:
                    branches,inbridges = self.expand_backward(cur)
                elif len(cur.out_nodes) == 1:
                    branches,inbridges = self.expand_forward(cur)
                else:
                    branches,inbridges = self.merge_join(cur)
                for branch in branches:
                    queue.appendleft(branch)
                for key, val in inbridges.items():
                    inbridge_list[key] = val
            elif len(sites) == 1:
                if self.next_is_stop(cur) or cur.is_bridge():
                    cur.split_node(sites[0], cleavage=True)
                    branches, inbridges = self.expand_forward(cur)
                    for branch in branches:
                        queue.appendleft(branch)
                    for key, val in inbridges.items():
                        inbridge_list[key] = val
                    continue
                if len(cur.out_nodes) == 1:
                    right = cur.split_node(sites[0], cleavage=True)
                    _,inbridges = self.expand_forward(cur)
                    queue.appendleft(list(right.out_nodes)[0])
                    for key, val in inbridges.items():
                        inbridge_list[key] = val
                    continue
                branches, inbridges = self.cross_join(cur, sites[0])
                for branch in branches:
                    queue.appendleft(branch)
                for key, val in inbridges.items():
                    inbridge_list[key] = val
            else:
                right = cur.split_node(sites[0], cleavage=True)
                if not cur.cleavage:
                    _,inbridges = self.expand_forward(cur)
                    for key, val in inbridges.items():
                        inbridge_list[key] = val
                site = right.seq.find_first_enzymatic_cleave_site(
                    rule=self.rule, exception=self.exception)
                right = right.split_node(site, cleavage=True)
                queue.appendleft(right)

    def update_orf(self, orf:List[int,int]) -> None:
        """ Update the orf list with the orf start position of the given node.
        """
        self.orfs.add(tuple(orf))

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
        cur = PVGCursor(None, self.root, True, [None,None], [])
        queue:Deque[Tuple[PVGNode,bool]] = deque([cur])
        peptide_pool = VariantPeptideDict(tx_id=self.id)
        traversal = PVGTraversal(
            check_variants=check_variants, check_orf=check_orf,
            miscleavage=miscleavage, queue=queue, pool=peptide_pool,
        )

        if self.has_known_orf():
            traversal.known_orf_aa = (
                int((self.known_orf[0] - self.known_reading_frame_index())/3),
                int((self.known_orf[1] - self.known_reading_frame_index())/3)
            )
            traversal.known_orf_tx = tuple(self.known_orf)

        while not traversal.is_done():
            cur = traversal.queue.pop()
            if cur.out_node is self.root and cur.out_node.seq is None:
                for out_node in cur.out_node.out_nodes:
                    cur = PVGCursor(cur.out_node, out_node, False,
                        cur.orf, [])
                    traversal.queue.appendleft(cur)
                continue

            if cur.out_node is self.stop:
                continue

            if self.has_known_orf():
                # add out nodes to the queue
                self.call_and_stage_known_orf(
                    cursor=cur,
                    traversal=traversal
                )
            else:
                self.call_and_stage_unknown_orf(
                    cursor=cur,
                    traversal=traversal
                )

        if check_orf:
            self.create_orf_id_map()

        return peptide_pool.get_peptide_sequences(
            keep_all_occurrence=keep_all_occurrence, orf_id_map=self.orf_id_map
        )

    def call_and_stage_known_orf(self, cursor:PVGCursor,
            traversal:PVGTraversal) -> None:
        """ For a given node in the graph, call miscleavage peptides if it
        is in CDS. For each of its outbond node, stage it until all inbond
        edges of the outbond node is visited. """
        target_node = cursor.out_node
        in_cds = cursor.in_cds
        orf = cursor.orf
        start_gain = cursor.start_gain
        if in_cds:
            node_copy = target_node.copy()

            stop_index = node_copy.seq.seq.find('*')

            if stop_index > -1:
                node_copy.truncate_right(stop_index)
                for out_node in copy.copy(node_copy.out_nodes):
                    node_copy.remove_out_edge(out_node)
            elif not node_copy.out_nodes:
                node_copy.truncated = True

            traversal.pool.add_miscleaved_sequences(
                node=node_copy,
                orf=tuple(orf),
                miscleavage=traversal.miscleavage,
                check_variants=traversal.check_variants,
                is_start_codon=False,
                additional_variants=start_gain
            )
            self.remove_node(node_copy)
            if stop_index > -1:
                in_cds = False
                orf = [None, None]

            for out_node in target_node.out_nodes:
                if out_node is self.stop:
                    continue
                if in_cds:
                    cur_start_gain = copy.copy(start_gain)
                    if not cur_start_gain:
                        new_start_gain = set()
                        for variant in out_node.variants:
                            if variant.variant.is_frameshifting():
                                new_start_gain.add(variant.variant)
                        for variant in target_node.variants:
                            if variant.variant.is_frameshifting():
                                new_start_gain.add(variant.variant)
                        cur_start_gain = list(new_start_gain)
                else:
                    cur_start_gain = []
                cur = PVGCursor(target_node, out_node, in_cds, orf,
                    cur_start_gain)
                traversal.stage(target_node, out_node, cur)
        else:
            if target_node.reading_frame_index != self.known_reading_frame_index():
                for out_node in target_node.out_nodes:
                    cur = PVGCursor(target_node, out_node, False, orf, [])
                    traversal.stage(target_node, out_node, cur)
                return

            if target_node.variants:
                start_index = target_node.seq.seq.find('M')
            else:
                start_index = target_node.seq.get_query_index(
                    traversal.known_orf_aa[0])
            if start_index == -1:
                for out_node in target_node.out_nodes:
                    cur = PVGCursor(target_node, out_node, False, orf, [])
                    traversal.stage(target_node, out_node, cur)
            else:
                node_copy = target_node.copy()
                stop_index = node_copy.seq.seq.find('*')
                if stop_index < start_index:
                    node_copy.truncate_left(start_index)
                    in_cds = True
                    more_stop_index = node_copy.seq.seq.find('*')
                    if more_stop_index > -1:
                        # checks if there are more stop after start.
                        node_copy.truncate_right(more_stop_index)
                        in_cds = False
                        for out_node in copy.copy(node_copy.out_nodes):
                            node_copy.remove_out_edge(out_node)
                else:
                    node_copy = node_copy[start_index:stop_index]
                    in_cds = False
                orf = traversal.known_orf_tx
                traversal.pool.add_miscleaved_sequences(
                    node=node_copy,
                    orf=tuple(orf),
                    miscleavage=traversal.miscleavage,
                    check_variants=traversal.check_variants,
                    is_start_codon=True,
                    additional_variants=[]
                )
                for out_node in target_node.out_nodes:
                    if out_node is not self.stop:
                        cur_start_gain = set()
                        for variant in out_node.variants:
                            if variant.variant.is_frameshifting():
                                cur_start_gain.add(variant.variant)
                        for variant in target_node.variants:
                            if variant.variant.is_frameshifting():
                                cur_start_gain.add(variant.variant)
                        cur_start_gain = list(cur_start_gain)
                        cur = PVGCursor(target_node, out_node, in_cds, orf,
                            cur_start_gain)
                        traversal.stage(target_node, out_node, cur)
                self.remove_node(node_copy)

    def call_and_stage_unknown_orf(self, cursor:PVGCursor,
            traversal:PVGTraversal) -> None:
        """ For a given node in the graph, call miscleavage peptides if it
        is in CDS. For each of its outbond node, stage it until all inbond
        edges of the outbond node is visited. """
        target_node = cursor.out_node
        in_cds = cursor.in_cds
        start_gain = cursor.start_gain
        orf = cursor.orf

        node_list:List[Tuple[PVGNode, List[int,int], bool]] = []
        trash = set()
        stop_index = target_node.seq.seq.find('*')
        start_index = target_node.seq.seq.find('M')

        if in_cds and stop_index == -1:
            node_list.append((target_node, orf, False))

        if in_cds and stop_index > -1:
            cur_copy = target_node.copy()
            cur_copy.truncate_right(stop_index)
            cur_copy.remove_out_edges()
            node_list.append((cur_copy, orf, False))
            trash.add(cur_copy)
            if stop_index > start_index:
                in_cds = False
                orf = [None, None]

        last_stop_index, last_start_index = -1, -1

        while stop_index > -1 or start_index > -1:
            if 0 < start_index < stop_index and (start_index !=
                    last_start_index or stop_index != last_stop_index):
                cur_copy = target_node[start_index:stop_index]
                orf_start = target_node.get_orf_start(start_index)
                cur_copy.remove_out_edges()
                cur_orf = [orf_start, None]
                self.add_stop(cur_copy)
                self.update_orf(cur_orf)
                node_list.append((cur_copy, cur_orf, True))
                trash.add(cur_copy)
            if start_index > -1 and stop_index == -1 \
                    and start_index != last_start_index:
                cur_copy = target_node.copy()
                cur_copy.truncate_left(start_index)
                orf_start = cur_copy.get_orf_start()
                cur_orf = [orf_start, None]
                self.update_orf(cur_orf)
                if not in_cds:
                    in_cds = True
                    orf = cur_orf
                node_list.append((cur_copy, cur_orf, True))
                trash.add(cur_copy)

            if start_index != -1:
                last_start_index = start_index
            if stop_index != -1:
                last_stop_index = stop_index
            # last_start_index, last_stop_index = start_index, stop_index

            if -1 < start_index < stop_index or (stop_index == -1 \
                    and start_index > -1):
                x = target_node.seq.seq[start_index+1:].find('M')
                start_index = x + start_index + 1 if x > -1 else x
            elif -1 < stop_index < start_index or (stop_index > -1 \
                    and start_index == -1):
                x = target_node.seq.seq[stop_index+1:].find('*')
                stop_index = x + stop_index + 1 if x > -1 else x

        for out_node in target_node.out_nodes:
            if out_node is not self.stop:
                out_node.orf = orf
                if target_node.is_bridge():
                    if not in_cds:
                        start_gain = []
                    elif last_start_index > -1:
                        # carry over variants from the target node to the next
                        # node if a start codon is found.
                        start_gain = target_node.get_variant_at(
                            start=last_start_index,
                            end=min(last_start_index + 3, len(target_node.seq.seq))
                        )
                    else:
                        start_gain = [v.variant for v in out_node.variants\
                            if v.variant.is_frameshifting()]
                else:
                    if last_stop_index > last_start_index:
                        start_gain = []

                cursor = PVGCursor(target_node, out_node, in_cds, orf,
                    start_gain)
                traversal.stage(target_node, out_node, cursor)

        for node, orf, is_start_codon in node_list:
            traversal.pool.add_miscleaved_sequences(
                node=node,
                orf=tuple(orf),
                miscleavage=traversal.miscleavage,
                check_variants=traversal.check_variants,
                is_start_codon=is_start_codon,
                additional_variants=start_gain
            )
        for node in trash:
            self.remove_node(node)

class PVGCursor():
    """ Helper class for cursors when graph traversal to call peptides. """
    def __init__(self, in_node:PVGNode, out_node:PVGNode, in_cds:bool,
            orf:List[int,int], start_gain:List[seqvar.VariantRecord]):
        """ constructor """
        self.in_node = in_node
        self.out_node = out_node
        self.in_cds = in_cds
        self.start_gain = start_gain
        self.orf = orf or [None, None]

class PVGTraversal():
    """ PVG Traversal. The purpose of this class is to facilitate the graph
    traversal to call variant peptides.
    """
    def __init__(self, check_variants:bool, check_orf:bool, miscleavage:int,
            pool:VariantPeptideDict, known_orf_tx:Tuple[int,int]=None,
            known_orf_aa:Tuple[int,int]=None,
            queue:Deque[PVGCursor]=None,
            stack:Dict[PVGNode, Dict[PVGNode, PVGCursor]]=None):
        """ constructor """
        self.check_variants = check_variants
        self.check_orf = check_orf
        self.miscleavage = miscleavage
        self.known_orf_tx = known_orf_tx or (None, None)
        self.known_orf_aa = known_orf_aa or (None, None)
        self.queue = queue or deque([])
        self.pool = pool
        self.stack = stack or {}

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

        if x.start_gain and not y.start_gain:
            return 1
        if not x.start_gain and y.start_gain:
            return -1
        if x.start_gain and y.start_gain:
            return -1 if sorted(x.start_gain)[0] > sorted(y.start_gain)[0] else 1
        return -1

    @staticmethod
    def cmp_known_orf_frame_shifted(x:PVGCursor, y:PVGCursor) -> bool:
        """ comparison for frameshifing nodes with known ORF """
        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1

        if x.start_gain and not y.start_gain:
            return -1
        if not x.start_gain and y.start_gain:
            return 1
        if x.start_gain and y.start_gain:
            return -1 if sorted(x.start_gain)[0] > sorted(y.start_gain)[0] else 1
        return -1

    @staticmethod
    def cmp_unknown_orf(x:PVGCursor, y:PVGCursor) -> bool:
        """ comparision for unkonwn ORF """
        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1

        if x.start_gain and not y.start_gain:
            return -1
        if not x.start_gain and y.start_gain:
            return 1
        if x.start_gain and y.start_gain:
            return -1 if sorted(x.start_gain)[0] > sorted(y.start_gain)[0] else 1
        return -1

    @staticmethod
    def cmp_unknown_orf_check_orf(x:PVGCursor, y:PVGCursor) -> bool:
        """ comparison for unknown ORF and check ORFs. """
        # pylint: disable=R0911
        if x.in_cds and not y.in_cds:
            return -1
        if not x.in_cds and y.in_cds:
            return 1

        if x.orf[0] is not None and (y.orf[0] is None or y.orf[0] == -1):
            return -1
        if (x.orf[0] is None or x.orf[0] == -1) and y.orf[0] is not None:
            return 1
        if x.orf[0] is not None and x.orf[0] != -1 and \
                y.orf[0] is not None and y.orf[0] != -1:
            if x.orf[0] > y.orf[0]:
                return -1
            if x.orf[0] < y.orf[0]:
                return 1

        if x.start_gain and not y.start_gain:
            return -1
        if not x.start_gain and y.start_gain:
            return 1
        if x.start_gain and y.start_gain:
            return -1 if sorted(x.start_gain)[0] > sorted(y.start_gain)[0] else 1
        return -1

    def stage(self, in_node:PVGNode, out_node:PVGNode,
            cursor:PVGCursor):
        """ When a node is visited through a particular edge during the variant
        peptide finding graph traversal, it is stage until all inbond edges
        are visited. """
        in_nodes = self.stack.setdefault(out_node, {})

        if in_node in in_nodes:
            return
        in_nodes[in_node] = cursor

        if len(in_nodes) != len(out_node.in_nodes):
            return

        curs = list(self.stack[out_node].values())
        if self.known_orf_aa[0] is not None:
            if out_node.reading_frame_index == self.known_reading_frame_index():
                func = self.cmp_known_orf_in_frame
            else:
                func = self.cmp_known_orf_frame_shifted
        elif self.check_orf:
            func = self.cmp_unknown_orf_check_orf
        else:
            func = self.cmp_unknown_orf

        curs.sort(key=cmp_to_key(func))

        self.queue.appendleft(curs[0])
