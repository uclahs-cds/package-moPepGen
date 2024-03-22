""" Module for peptide variation graph """
from __future__ import annotations
import copy
from typing import Callable, FrozenSet, Iterable, Set, Deque, Dict, List, Tuple, TYPE_CHECKING
from collections import deque
from functools import cmp_to_key
from Bio.Seq import Seq
from moPepGen import aa, params
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.svgraph.SubgraphTree import SubgraphTree
from moPepGen.svgraph.VariantPeptideDict import VariantPeptideDict
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.svgraph.PVGOrf import PVGOrf
from moPepGen.svgraph.PVGNodeCollapser import PVGNodeCollapser, PVGNodePopCollapser


if TYPE_CHECKING:
    from .VariantPeptideDict import AnnotatedPeptideLabel
    from moPepGen.circ import CircRNAModel
    from moPepGen.params import CleavageParams

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
            known_orf:List[int,int], cleavage_params:CleavageParams=None,
            orfs:Set[Tuple[int,int]]=None, reading_frames:List[PVGNode]=None,
            orf_id_map:Dict[int,str]=None, cds_start_nf:bool=False,
            hypermutated_region_warned:bool=False, denylist:Set[str]=None,
            global_variant:VariantRecord=None, subgraphs:SubgraphTree=None,
            gene_id:str=None):
        """ Construct a PeptideVariantGraph """
        self.root = root
        self.id = _id
        self.known_orf = known_orf
        stop_seq = aa.AminoAcidSeqRecordWithCoordinates(
            seq=Seq('*'), locations=[]
        )
        self.stop = PVGNode(stop_seq, None, subgraph_id=self.id)
        self.orfs = orfs or set()
        self.reading_frames = reading_frames or [None, None, None]
        self.orf_id_map = orf_id_map or {}
        self.cds_start_nf = cds_start_nf
        self.hypermutated_region_warned = hypermutated_region_warned
        self.denylist = denylist or set()
        self.cleavage_params = cleavage_params or params.CleavageParams()
        self.global_variant = global_variant
        self.subgraphs = subgraphs
        self.gene_id = gene_id

    def add_stop(self, node:PVGNode):
        """ Add the stop node after the specified node. """
        node.add_out_edge(self.stop)

    def next_is_stop(self, node:PVGNode) -> bool:
        """ Checks if stop is linked as a outbound node. """
        return len(node.out_nodes) == 1 and self.stop in node.out_nodes

    def is_circ_rna(self) -> bool:
        """ Checks if this is a circRNA """
        return self.global_variant and self.global_variant.is_circ_rna()

    def find_node(self, func:Callable) -> List[PVGNode]:
        """ Find node with given sequence. """
        queue:Deque[PVGNode] = deque([self.root])
        visited:Set[PVGNode] = set()
        targets:List[PVGNode] = []
        while queue:
            cur = queue.popleft()
            if cur in visited:
                continue
            visited.add(cur)
            # pylint: disable=W0703
            try:
                if func(cur):
                    targets.append(cur)
            except Exception:
                pass
            finally:
                for out_node in cur.out_nodes:
                    if not out_node in visited:
                        queue.append(out_node)
        return targets

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

    def remove_edge(self, in_node:PVGNode, out_node:PVGNode) -> None:
        """ Remove an edge """
        in_node.remove_out_edge(out_node)
        if not in_node.out_nodes:
            self.remove_node(in_node)
        if not out_node.in_nodes:
            self.remove_node(out_node)

    def known_reading_frame_index(self) -> int:
        """ Get the readign frame index of the known ORF. """
        return self.known_orf[0] % 3

    def has_known_orf(self) -> bool:
        """ Check if it has a known ORF """
        return self.known_orf[0] is not None

    def node_is_subgraph_end(self, node:PVGNode) -> bool:
        """ check if a node is the last node of a subgraph """
        return node.level > self.root.level and \
            any(x.subgraph_id != node.subgraph_id and x.level < node.level
                and x is not self.stop for x in node.out_nodes)

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
        exception_sites = node.seq.get_enzymatic_cleave_exception_sites(
            self.cleavage_params.exception
        )
        sites = node.seq.find_all_cleave_and_stop_sites_with_range(
            rule=self.cleavage_params.enzyme,
            exception=self.cleavage_params.exception,
            exception_sites=exception_sites
        )
        sites = deque(sites)
        shift = 0
        while len(sites) > 0:
            s, r = sites.popleft()
            shifted_site = s - shift
            if r:
                r = (r[0] - shift, r[1] - shift)
            shift = s
            node = node.split_node(shifted_site, cleavage=True, cleavage_range=r)
        return first_node if return_first else node

    def count_nodes(self):
        """ Count nodes """
        queue = deque([self.root])
        visited = set()
        k = 0

        while queue:
            cur = queue.pop()
            if cur in visited:
                continue
            k += 1
            visited.add(cur)
            for node in cur.out_nodes:
                if node in visited:
                    continue
                if node is None:
                    continue
                queue.appendleft(node)
        return k

    def find_nodes_with_seq(self, seq:str) -> List[PVGNode]:
        """ find all nodes with the given sequence """
        queue = deque(self.root.out_nodes)
        visited = set([self.root])
        found = []
        while queue:
            cur = queue.popleft()
            if cur in visited:
                continue
            if cur.seq.seq == seq:
                found.append(cur)
            visited.add(cur)
            for out_node in cur.out_nodes:
                queue.append(out_node)
        return found

    def nodes_have_too_many_variants(self, nodes:Iterable[PVGNode]) -> bool:
        """ Check the total number of variants of given nodes """
        if self.cleavage_params.max_variants_per_node == -1:
            return False
        variants = set()
        for node in nodes:
            for variant in node.variants:
                variants.add(variant.variant)
        return len(variants) > self.cleavage_params.max_variants_per_node

    def find_routes_for_merging(self, node:PVGNode, cleavage:bool=False
            ) -> Tuple[Set[Tuple[PVGNode]], Set[PVGNode]]:
        """ Find all start and end nodes for merging.

        Args:
            node (PVGNode)
            cleavage (bool): When True, the input node is treated as the start
                node, so it's in-bond nodes are not looked.

        Returns:
            A tuple of two. The first is lists of nodes that should be merged,
            and the second is all the original nodes that should be removed.
        """
        routes:Set[Tuple[PVGNode]] = set()
        visited:Set[PVGNode] = {node}
        for out_node in node.out_nodes:
            if out_node.is_bridge() or self.node_is_subgraph_end(out_node):
                continue
            is_sharing_downstream = any(any(y is not out_node and y in node.out_nodes \
                for y in x.in_nodes) for x in out_node.out_nodes if x is not self.stop)
            if is_sharing_downstream:
                continue
            s, r = out_node.seq.find_first_cleave_or_stop_site_with_range(
                rule=self.cleavage_params.enzyme,
                exception=self.cleavage_params.exception
            )
            if s > -1:
                out_node.split_node(s, cleavage=True, pre_cleave=True, cleavage_range=r)

        if cleavage:
            for out_node in node.out_nodes:
                count = 0
                route = (node, out_node)
                for out_node_2 in out_node.out_nodes:
                    if node.is_inbond_of(out_node_2):
                        count += 1
                        new_route = (*route, out_node_2)
                        routes.add(new_route)
                if count == 0:
                    routes.add(route)
                    visited.add(node)
                elif count == len(out_node.out_nodes):
                    visited.add(node)

                for in_node in out_node.in_nodes:
                    if in_node is not node and not node.is_inbond_of(in_node):
                        route = (in_node, out_node)
                        routes.add(route)
                        visited.add(in_node)
        else:
            for out_node in node.out_nodes:
                visited.add(out_node)
                for in_node in node.in_nodes:
                    route = (in_node, node, out_node)
                    routes.add(route)
                    visited.add(in_node)

                for in_node in out_node.in_nodes:
                    if in_node is not node:
                        route = (in_node, out_node)
                        routes.add(route)
                        visited.add(in_node)

        return routes, visited

    def merge_nodes_routes(self, routes:Set[Tuple[PVGNode]],
            reading_frame_index:int
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

        for route in routes:
            for i,node in enumerate(route):
                node_is_bridge = node.is_bridge()
                if i == 0:
                    new_node = node.copy(in_nodes=False, out_nodes=False)
                else:
                    new_node.append_right(node)
                    if new_node.level < node.level:
                        new_node.subgraph_id = node.subgraph_id
                        new_node.level = node.level
                    if node_is_bridge:
                        new_node.was_bridge = True
                    trash.add((route[i-1], node))

            for in_node in route[0].in_nodes:
                in_node.add_out_edge(new_node)

            for out_node in route[-1].out_nodes:
                new_node.add_out_edge(out_node)

            new_nodes.add(new_node)

            is_in_bridge = route[0].reading_frame_index != reading_frame_index
            is_out_bridge = route[-1].is_bridge() and \
                route[-1].reading_frame_index == reading_frame_index
            if is_in_bridge and not is_out_bridge or self.node_is_subgraph_end(route[0]):
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
        downstream_2_nodes:Dict[FrozenSet[PVGNode], List[PVGNode]] = {}
        for node in nodes:
            out_nodes = frozenset(node.out_nodes)
            if out_nodes in downstream_2_nodes:
                downstream_2_nodes[out_nodes].append(node)
            else:
                downstream_2_nodes[out_nodes] = [node]
        for node in nodes:
            if node.seq.seq == '*' and not node.out_nodes:
                continue
            if node.get_last_rf_index() != reading_frame_index \
                    and len(node.get_out_nodes()) == 1 \
                    and not node.has_exclusive_outbond_node() \
                    and not all(x in nodes for x in node.get_out_nodes()[0].get_in_nodes()) \
                    and not len(node.get_out_nodes()[0].get_out_nodes()) == 0 \
                    and not node.get_out_nodes()[0].get_out_nodes()[0].seq.seq == '*':
                continue
            is_deletion_only_end = any(x.variant.type == 'Deletion' for x in node.variants) \
                and len(node.out_nodes) == 1 \
                and len(node.get_out_nodes()[0].in_nodes) > 1 \
                and not all(any(v.variant.type == 'Deletion' for v in x.variants)
                    for x in node.get_out_nodes()[0].get_in_nodes())
            if node.is_bridge() \
                    or self.node_is_subgraph_end(node) \
                    or is_deletion_only_end:
                continue

            is_sharing_downstream = any(any(y is not node and y in nodes \
                for y in x.in_nodes) for x in node.out_nodes if x is not self.stop)

            if not is_sharing_downstream:
                is_single_out = len(node.out_nodes) == 1 and \
                    len(list(node.out_nodes)[0].in_nodes) == 1 and \
                    node.cleavage and list(node.out_nodes)[0].cleavage
                if is_single_out:
                    out_node = list(node.out_nodes)[0]
                    if out_node.is_bridge() or self.node_is_subgraph_end(out_node):
                        continue
                    downstreams.add(out_node)
                else:
                    downstreams.add(node)
            else:
                if len(node.out_nodes) > 1:
                    downstreams.update(downstream_2_nodes[frozenset(node.out_nodes)])
                    continue
                for downstream in node.out_nodes:
                    if downstream is self.stop:
                        continue
                    downstreams.add(downstream)
        return downstreams

    def update_unique_nodes(self, node:PVGNode, unique_nodes:Set[PVGNode]):
        """ For a given node and a set of unique nodes, if the node has an
        identical mate in the set, collpase them into one and keep all
        inbond node to it. Otherwise, add the node to the set. """
        node_collapsed  = False
        for unique_node in copy.copy(unique_nodes):
            if unique_node.is_identical(node):
                node_collapsed = True
                if node.is_less_mutated_than(unique_node):
                    unique_nodes.remove(unique_node)
                    unique_nodes.add(node)
                    unique_node.transfer_in_nodes_to(node)
                    self.remove_node(unique_node)
                else:
                    node.transfer_in_nodes_to(unique_node)
                    self.remove_node(node)
                break
        if not node_collapsed:
            unique_nodes.add(node)

    def collapse_in_nodes(self, node:PVGNode) -> Set[PVGNode]:
        """ Collapse all inbond nodes if they are identical. """
        unique_nodes:Set[PVGNode] = set()
        for in_node in copy.copy(node.in_nodes):
            self.update_unique_nodes(in_node, unique_nodes)
        return unique_nodes

    def collapse_end_nodes(self, nodes:Iterable[PVGNode]):
        """ Collapse nodes inthey are identifical. This function is called
        after 'routes' are merged. Then redundant end nodes are collapsed. """
        group:Dict[Tuple[PVGNode],PVGNodeCollapser] = {}
        for node in nodes:
            if len(node.out_nodes) == 1 and self.next_is_stop(node):
                continue
            out_nodes = tuple(node.out_nodes)
            collapser = group.setdefault(out_nodes, PVGNodeCollapser())
            redundant_node = collapser.collapse(node)
            if redundant_node:
                self.remove_node(redundant_node)
        collapsed_nodes = set(nodes)
        for node in copy.copy(collapsed_nodes):
            if node.is_orphan():
                collapsed_nodes.remove(node)
        return collapsed_nodes

    def collapse_nodes_backward(self, nodes:Iterable[PVGNode], heads:Set[PVGNode]
            ) -> Set[PVGNode]:
        """ Collapse equivalent nodes in a variant bubble from end to start.
        Equivalent nodes are those with same sequence and share the same
        outgoing nodes. """
        final_nodes = self.collapse_end_nodes(nodes)
        queue = deque(final_nodes)
        while queue:
            cur = queue.popleft()
            if any(x in heads for x in cur.get_in_nodes()):
                continue
            collapsed_nodes = self.collapse_end_nodes(cur.get_in_nodes())
            for collapsed_node in collapsed_nodes:
                queue.append(collapsed_node)
        return final_nodes

    def pop_collapse_end_nodes(self, nodes:Iterable[PVGNode]):
        """ This function aims at resolving the issue that too may nodes are
        generated when making the cleavage graph. For nodes that share the
        outbond nodes, the last X number of amino acids are popped as separate
        nodes, and collapsed if they have the same sequence, and the least
        variant node is kept. """
        group:Dict[FrozenSet[PVGNode], PVGNodePopCollapser] = {}
        collapsed_nodes:List[PVGNode] = []
        for node in nodes:
            node_len = len(node.seq.seq)
            if node_len <= self.cleavage_params.naa_to_collapse:
                collapsed_nodes.append(node)
                continue
            new_node = node.split_node(
                index=node_len - self.cleavage_params.naa_to_collapse,
                cleavage=True,
                pop_collapse=True
            )
            collapsed_nodes.append(new_node)
            out_nodes = frozenset(new_node.out_nodes)
            collapser = group.setdefault(out_nodes, PVGNodePopCollapser())
            redundant_node = collapser.collapse(new_node)
            if redundant_node:
                self.remove_node(redundant_node)
        collapsed_nodes = set(collapsed_nodes)
        for node in copy.copy(collapsed_nodes):
            if node.is_orphan():
                collapsed_nodes.remove(node)
        return collapsed_nodes

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
        routes, _ = self.find_routes_for_merging(node, True)

        heads = set()
        for x in routes:
            heads.update(x[0].get_in_nodes())

        new_nodes, inbridge_list = self.merge_nodes_routes(routes, reading_frame_index)
        new_nodes = self.collapse_nodes_backward(new_nodes, heads)
        if len(new_nodes) > self.cleavage_params.min_nodes_to_collapse:
            new_nodes = self.pop_collapse_end_nodes(new_nodes)
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
        routes:Set[PVGNode] = set()
        trash:Set[PVGNode] = set([node])
        reading_frame_index = node.reading_frame_index
        for in_node in node.in_nodes:
            routes.add((in_node, node))
            trash.add(in_node)

        heads = set()
        for x in routes:
            heads.update(x[0].get_in_nodes())

        new_nodes, inbridge_list = self.merge_nodes_routes(routes, reading_frame_index)
        new_nodes = self.collapse_nodes_backward(new_nodes, heads)
        # new_nodes = self.pop_collapse_end_nodes(new_nodes)
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
        routes, _ = self.find_routes_for_merging(node)

        heads = set()
        for x in routes:
            heads.update(x[0].get_in_nodes())

        new_nodes, inbridge_list = self.merge_nodes_routes(routes, reading_frame_index)
        new_nodes = self.collapse_nodes_backward(new_nodes, heads)
        if len(new_nodes) > self.cleavage_params.min_nodes_to_collapse:
            new_nodes = self.pop_collapse_end_nodes(new_nodes)
        downstreams = self.move_downstreams(new_nodes, reading_frame_index)
        return downstreams, inbridge_list

    def cross_join(self, node:PVGNode, site:int, cleavage_range:Tuple[int, int]=None) -> T:
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
        head = node.split_node(site, cleavage=True, cleavage_range=cleavage_range)

        # there shouldn't have any bridge out in inbond nodes
        left_nodes, _ = self.expand_forward(node)
        self.collapse_end_nodes(left_nodes)
        # self.pop_collapse_end_nodes(new_nodes)

        routes, _ = self.find_routes_for_merging(head, True)

        heads = set()
        for x in routes:
            heads.update(x[0].get_in_nodes())

        new_nodes, inbridge_list = self.merge_nodes_routes(routes, reading_frame_index)
        new_nodes = self.collapse_nodes_backward(new_nodes, heads)
        if len(new_nodes) > self.cleavage_params.min_nodes_to_collapse:
            new_nodes = self.pop_collapse_end_nodes(new_nodes)
        downstreams = self.move_downstreams(new_nodes, reading_frame_index)
        return downstreams, inbridge_list

    def create_cleavage_graph(self) -> None:
        """ Form a cleavage graph from a variant graph. After calling this
        method, each in the graph should represent a cleavage in the
        sequence of the pull peptide of reference and variant.

        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
        """
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

            first_site = cur.seq.find_first_cleave_or_stop_site(
                rule=self.cleavage_params.enzyme,
                exception=self.cleavage_params.exception
            )

            if cur.is_already_cleaved() and first_site == -1:
                continue

            if len(cur.in_nodes) == 1:
                downstreams, inbridges = self.fit_into_cleavages_single_upstream(cur)
            else:
                downstreams, inbridges = \
                    self.fit_into_cleavage_multiple_upstream(cur)

            for downstream in downstreams:
                queue.appendleft(downstream)
            for key, val in inbridges.items():
                inbridge_list[key] = val

    def fit_into_cleavages_single_upstream(self, cur:PVGNode) -> T:
        """ Fit node into cleavage sites when it has a single upstream node """
        cur = self.cleave_if_possible(cur)
        if self.next_is_stop(cur):
            if len(cur.out_nodes) > 1:
                return self.expand_backward(cur)
            return set(), {}

        if cur.is_bridge() or self.node_is_subgraph_end(cur):
            return set(), {}

        i = cur.reading_frame_index

        if len(cur.out_nodes) == 1:
            downstream = list(cur.out_nodes)[0]
            s, r = downstream.seq.find_first_cleave_or_stop_site_with_range(
                rule=self.cleavage_params.enzyme,
                exception=self.cleavage_params.exception
            )
            if s > -1:
                downstream.split_node(s, True, cleavage_range=r)
            # It's important that if the only downstream is bridge, we only
            # do a simple expand forward. Issue #
            if len(downstream.out_nodes) == 1 or downstream.is_bridge():
                branches, inbridges = self.expand_forward(downstream)
            else:
                branches, inbridges = self.merge_join(downstream)
            branches = {x for x in branches if x.reading_frame_index == i}
        else:
            branches, inbridges = self.expand_backward(cur)
            branches = {x for x in branches if x.reading_frame_index == i and\
                not (x.cleavage and (x.is_bridge() or self.node_is_subgraph_end(x)))}
        return branches, inbridges

    def fit_into_cleavage_multiple_upstream(self, cur:PVGNode) -> T:
        """ Fit a node into cleavage sites when it has multiple upstream nodes
        """
        branches:Set[PVGNode] = set()
        inbridges:Dict[PVGNode,List[PVGNode]] = {}

        sites = cur.seq.find_all_cleave_and_stop_sites_with_range(
            rule=self.cleavage_params.enzyme,
            exception=self.cleavage_params.exception
        )

        if len(sites) == 0:
            if self.next_is_stop(cur) or cur.is_bridge() or \
                    self.node_is_subgraph_end(cur):
                if not cur.cleavage:
                    _,inbridges = self.expand_forward(cur)
                return branches, inbridges
            if cur.cleavage:
                branches,inbridges = self.expand_backward(cur)
            elif len(cur.out_nodes) == 1:
                branches,inbridges = self.merge_join(cur)
            else:
                branches,inbridges = self.merge_join(cur)

        elif len(sites) == 1:
            s, r = sites[0]
            if self.next_is_stop(cur) or cur.is_bridge() or \
                    self.node_is_subgraph_end(cur):
                cur.split_node(s, cleavage=True, cleavage_range=r)
                branches, inbridges = self.expand_forward(cur)
            elif len(cur.out_nodes) == 1:
                right = cur.split_node(s, cleavage=True, cleavage_range=r)
                _,inbridges = self.expand_forward(cur)
                branches = {right}
            else:
                branches, inbridges = self.cross_join(cur, s, cleavage_range=r)

        else:
            if self.next_is_stop(cur) or cur.is_bridge() or \
                    self.node_is_subgraph_end(cur):
                cur = self.cleave_if_possible(cur, return_first=True)
                if not cur.cleavage:
                    _,inbridges = self.expand_forward(cur)
            else:
                s, r = sites[0]
                right = cur.split_node(s, cleavage=True, cleavage_range=r)
                if not cur.cleavage:
                    _,inbridges = self.expand_forward(cur)
                s, r = right.seq.find_first_cleave_or_stop_site_with_range(
                    rule=self.cleavage_params.enzyme,
                    exception=self.cleavage_params.exception
                )
                right = right.split_node(s, cleavage=True, cleavage_range=r)
                branches = {right}

        return branches, inbridges

    def update_orf(self, orf:List[int,int]) -> None:
        """ Update the orf list with the orf start position of the given node.
        """
        self.orfs.add(tuple(orf))

    def create_orf_id_map(self) -> None:
        """ Creates a map for ORF start site and its ID (ORF1, ORF2, etc) """
        orfs = list(self.orfs)
        orfs.sort()
        self.orf_id_map = {v:f"ORF{i+1}" for i,v in enumerate(orfs)}

    def call_variant_peptides(self, check_variants:bool=True,
            check_orf:bool=False, keep_all_occurrence:bool=True, denylist:Set[str]=None,
            circ_rna:CircRNAModel=None, orf_assignment:str='max',
            backsplicing_only:bool=False, truncate_sec:bool=False, w2f:bool=False,
            check_external_variants:bool=True, find_ass:bool=False
            ) -> Dict[Seq, List[AnnotatedPeptideLabel]]:
        """ Walk through the graph and find all noncanonical peptides.

        Args:
        - miscleavage (int): Number of miscleavages allowed.
        - check_variants (bool): When true, only peptides that carry at
            least 1 variant are kept. When false, all unique peptides
            are reported.
        - check_orf (bool): When true, the ORF ID will be added to the
            variant peptide label.
        - keep_all_occurrence (bool): Whether to keep all occurences of the
            peptide within the graph/transcript.
        - denylist (Set[str]): Denylist of peptides to be called as variant
            peptides.
        - circ_rna (CircRNAModel): The circRNA from which variant peptides are
            called.
        - orf_assignment (str): ORF assignment strategy. Options: ['max', 'min']
        - backsplicing_only (bool): Whether to only output variant peptides
            spanning the backsplicing site.
        - truncate_sec (bool): Whether to consider selenocysteine termination
            events.
        - w2f (bool): Whether to consider W>F (tryptophan to phenylalanine)
            substituants.
        - check_external_variants (bool): When set to `False`, peptides not
            harbouring any external variants will also be outputted. This can
            be used when calling noncanonical peptides from noncoding
            transcripts.
        - find_ass (bool): Whether to find alternative start sites.
        """
        if self.is_circ_rna() and not circ_rna:
            raise ValueError('`circ_rna` must be given.')

        self.denylist = denylist or set()
        cur = PVGCursor(None, self.root, True, [], [])
        queue:Deque[Tuple[PVGNode,bool]] = deque([cur])
        peptide_pool = VariantPeptideDict(
            tx_id=self.id,
            global_variant=self.global_variant,
            gene_id=self.gene_id,
            truncate_sec=truncate_sec,
            w2f=w2f,
            check_external_variants=check_external_variants,
            cleavage_params=self.cleavage_params,
            check_orf=check_orf
        )
        traversal = PVGTraversal(
            check_variants=check_external_variants,
            check_orf=check_orf, queue=queue, pool=peptide_pool,
            circ_rna=circ_rna, orf_assignment=orf_assignment,
            backsplicing_only=backsplicing_only,
            find_ass=find_ass
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
                        cur.orfs, [])
                    traversal.queue.appendleft(cur)
                continue

            if cur.out_node is self.stop:
                continue

            if self.is_circ_rna() and cur.out_node.is_hybrid_node(self.subgraphs):
                self.call_and_stage_silently(
                    cursor=cur, traversal=traversal
                )
            elif self.has_known_orf() and not traversal.find_ass:
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

        peptide_pool.translational_modification(w2f, self.denylist)

        return peptide_pool.get_peptide_sequences(
            keep_all_occurrence=keep_all_occurrence,
            orf_id_map=self.orf_id_map,
            check_variants=check_variants
        )

    def call_and_stage_known_orf(self, cursor:PVGCursor,
            traversal:PVGTraversal) -> None:
        """ For a given node in the graph, call miscleavage peptides if it
        is in CDS. For each of its outbond node, stage it until all inbond
        edges of the outbond node is visited. """
        if cursor.in_cds:
            self.call_and_stage_known_orf_in_cds(
                cursor=cursor, traversal=traversal
            )
        else:
            self.call_and_stage_known_orf_not_in_cds(
                cursor=cursor, traversal=traversal
            )

    def call_and_stage_known_orf_in_cds(self, cursor:PVGCursor,
            traversal:PVGTraversal) -> None:
        """ Known ORF in CDS branch """
        target_node = cursor.out_node
        in_cds = cursor.in_cds
        orfs = [orf.copy() for orf in cursor.orfs]

        is_stop = target_node.seq.seq == '*'

        if is_stop:
            in_cds = False
            orfs = []
        elif not target_node.npop_collapsed:
            node_copy = target_node.copy(in_nodes=False)

            additional_variants = cursor.cleavage_gain
            upstream_indels = target_node.upstream_indel_map.get(cursor.in_node)
            if upstream_indels:
                additional_variants += upstream_indels

            traversal.pool.add_miscleaved_sequences(
                node=node_copy,
                orfs=orfs,
                cleavage_params=self.cleavage_params,
                check_variants=traversal.check_variants,
                is_start_codon=False,
                additional_variants=additional_variants,
                denylist=self.denylist,
                leading_node=target_node,
                subgraphs=self.subgraphs,
                backsplicing_only=traversal.backsplicing_only
            )
            self.remove_node(node_copy)

        cleavage_gain = target_node.get_cleavage_gain_variants()

        for out_node in target_node.out_nodes:
            if out_node is self.stop:
                continue
            cur_orfs = [orf.copy() for orf in orfs]
            if in_cds:
                cur_start_gain = copy.copy(cur_orfs[0].start_gain)
                if not cur_start_gain:
                    cur_start_gain = set()

                    for variant in out_node.variants:
                        if variant.variant.is_frameshifting():
                            cur_start_gain.add(variant.variant)

                    for variant in target_node.variants:
                        if variant.variant.is_frameshifting():
                            cur_start_gain.add(variant.variant)

                    upstream_indels = target_node.upstream_indel_map.get(cursor.in_node)
                    if upstream_indels:
                        for variant in upstream_indels:
                            if variant.is_frameshifting():
                                cur_start_gain.add(variant)

                    stop_index = self.known_orf[1]
                    stop_lost = target_node.get_stop_lost_variants(stop_index)
                    cur_start_gain.update(stop_lost)
                cur_cleavage_gain = copy.copy(cleavage_gain)
                cleavage_gain_down = out_node.get_cleavage_gain_from_downstream()
                cur_cleavage_gain.extend(cleavage_gain_down)
                cur_orfs[0].start_gain = cur_start_gain
            else:
                cur_cleavage_gain = None
            cur = PVGCursor(target_node, out_node, in_cds, cur_orfs, cur_cleavage_gain)
            traversal.stage(target_node, out_node, cur)

    def call_and_stage_known_orf_not_in_cds(self, cursor:PVGCursor,
            traversal:PVGTraversal) -> None:
        """ Kown ORF not in CDS branch """
        target_node = cursor.out_node
        in_cds = cursor.in_cds
        orfs = []
        start_gain = set()
        if target_node.reading_frame_index != self.known_reading_frame_index():
            for out_node in target_node.out_nodes:
                cur = PVGCursor(target_node, out_node, False, orfs)
                traversal.stage(target_node, out_node, cur)
            return

        start_index = target_node.seq.get_query_index(
            ref_index=traversal.known_orf_aa[0],
            seqname=self.id,
            reading_frame=traversal.known_orf_tx[0] % 3
        )
        if start_index == -1:
            for out_node in target_node.out_nodes:
                cur = PVGCursor(
                    in_node=target_node, out_node=out_node,
                    in_cds=False, orfs=orfs
                )
                traversal.stage(target_node, out_node, cur)
        else:
            start_gain.update(target_node.get_variants_at(start_index))
            additional_variants = []
            node_copy = target_node.copy(in_nodes=False)
            in_cds = True
            node_copy.truncate_left(start_index)
            orf = PVGOrf(orf=list(traversal.known_orf_tx), start_gain=start_gain,
                start_node=target_node, subgraph_id=self.id, node_offset=start_index)
            traversal.pool.add_miscleaved_sequences(
                node=node_copy,
                orfs=[orf],
                cleavage_params=self.cleavage_params,
                check_variants=traversal.check_variants,
                is_start_codon=True,
                additional_variants=additional_variants,
                denylist=self.denylist,
                leading_node=target_node,
                subgraphs=self.subgraphs,
                backsplicing_only=traversal.backsplicing_only
            )
            cleavage_gain = target_node.get_cleavage_gain_variants()
            for out_node in target_node.out_nodes:
                if out_node is not self.stop:
                    cur_start_gain = set()
                    for variant in out_node.variants:
                        if variant.variant.is_frameshifting():
                            cur_start_gain.add(variant.variant)
                    for variant in target_node.variants:
                        if variant.variant.is_frameshifting():
                            cur_start_gain.add(variant.variant)
                        if variant.location.start > start_index \
                                and variant.is_stop_altering:
                            cur_start_gain.add(variant.variant)
                    cur_cleavage_gain = copy.copy(cleavage_gain)
                    cleavage_gain_down = out_node.get_cleavage_gain_from_downstream()
                    cur_cleavage_gain.extend(cleavage_gain_down)
                    cur_orf = orf.copy()
                    cur_orf.start_gain = cur_start_gain
                    cur = PVGCursor(
                        in_node=target_node, out_node=out_node, in_cds=in_cds,
                        orfs=[cur_orf],  cleavage_gain=cur_cleavage_gain
                    )
                    traversal.stage(target_node, out_node, cur)
            self.remove_node(node_copy)

    @staticmethod
    def call_and_stage_silently(cursor:PVGCursor, traversal:PVGTraversal):
        """ This is called when the cursor node is invalid. """
        target_node = cursor.out_node
        finding_start_site = cursor.finding_start_site

        for node in target_node.out_nodes:
            cursor = PVGCursor(target_node, node, False, [],
                [], finding_start_site)
            traversal.stage(target_node, node, cursor)

    def call_and_stage_unknown_orf(self, cursor:PVGCursor,
            traversal:PVGTraversal) -> None:
        """ For a given node in the graph, call miscleavage peptides if it
        is in CDS. For each of its outbond node, stage it until all inbond
        edges of the outbond node is visited. """
        target_node = cursor.out_node
        in_cds = cursor.in_cds
        orfs = [orf.copy() for orf in cursor.orfs]
        finding_start_site = cursor.finding_start_site

        node_list:List[Tuple[PVGNode, List[int,int], bool, List[VariantRecord]]] = []
        trash = set()
        is_stop = target_node.seq.seq == '*'

        cleavage_gain_down = target_node.get_cleavage_gain_from_downstream()

        if is_stop:
            in_cds = False
            orfs = []

        if in_cds and not target_node.npop_collapsed:
            cur_copy = target_node.copy(in_nodes=False)
            cur_orfs = copy.copy(orfs)
            if self.is_circ_rna():
                additional_variants = []
            else:
                additional_variants = cursor.cleavage_gain

            node_list.append((cur_copy, cur_orfs, False, additional_variants))
            trash.add(cur_copy)

        # if the current node contains the actual fusion variant, stop looking
        # for further start sites.
        if finding_start_site:
            for variant in target_node.variants:
                if variant.variant.is_real_fusion and variant.not_cleavage_altering():
                    finding_start_site = False
                    real_fusion_position = variant.location.start
                    break

        start_indices = []
        if cursor.finding_start_site:
            start_indices = target_node.seq.find_all_start_sites()
            if not finding_start_site:
                start_indices = [x for x in start_indices if x <= real_fusion_position]
            for start_index in start_indices:
                cur_copy = target_node.copy(in_nodes=False)
                cur_copy.truncate_left(start_index)
                orf_start = cur_copy.get_orf_start()
                orf_subgraph_id = cur_copy.get_subgraph_id_at(0)
                cur_orf = [orf_start, None]
                self.update_orf(cur_orf)
                in_cds = True
                orf = PVGOrf(
                    orf=cur_orf, start_gain=set(), start_node=target_node,
                    subgraph_id=orf_subgraph_id, node_offset = start_index
                )
                orfs.append(orf)
                if self.is_circ_rna():
                    additional_variants = []
                    orf.cleavage_gain = set(cleavage_gain_down)
                else:
                    additional_variants = copy.copy(cleavage_gain_down)
                cur_orfs = [orfs[0]] if traversal.orf_assignment == 'min' else [orfs[-1]]
                node_list.append((cur_copy, cur_orfs, True, additional_variants))
                trash.add(cur_copy)

        cleavage_gain = target_node.get_cleavage_gain_variants()

        for out_node in target_node.out_nodes:
            cur_in_cds = in_cds
            if out_node is self.stop:
                continue

            if not in_cds:
                start_gain = set()
                cur_orfs = []
            elif start_indices:
                cur_orf = orf.copy()
                # carry over variants from the target node to the next
                # node if a start codon is found.
                start_gain = target_node.get_variants_at(
                    start=start_indices[-1],
                    end=min(start_indices[-1] + 1, len(target_node.seq.seq)),
                    upstream_cleavage_altering=False,
                    downstream_cleavage_altering=False
                )
                fs_variants = target_node.get_variants_at(
                    start=start_indices[-1], end=-1,
                    upstream_cleavage_altering=False,
                    downstream_cleavage_altering=False
                )
                start_gain += [x for x in fs_variants if x.is_frameshifting()]
                start_gain = [x for x in start_gain if not x.is_circ_rna()]
                cur_orf.start_gain = set(start_gain)
                cur_orfs = [cur_orf]
            else:
                cur_orfs = [orf.copy() for orf in orfs]
                for orf in cur_orfs:
                    if self.is_circ_rna():
                        start_gain = [v.variant for v in target_node.variants
                            if not v.variant.is_circ_rna() and v.not_cleavage_altering()]
                        orf.start_gain.update(start_gain)
                    elif not orf.start_gain:
                        start_gain = [v.variant for v in target_node.variants
                            if v.variant.is_frameshifting()
                                and v.not_cleavage_altering()
                                and not v.variant.is_circ_rna()]
                        orf.start_gain.update(start_gain)

            # Add stop altering mutations
            for variant in target_node.variants:
                if variant.is_stop_altering and variant.not_cleavage_altering():
                    if start_indices and variant.location.end - 1 < start_indices[-1]:
                        continue
                    for x in cur_orfs:
                        x.start_gain.add(variant.variant)

            if self.is_circ_rna():
                cur_cleavage_gain = []
            else:
                cur_cleavage_gain = copy.copy(cleavage_gain)

            if self.is_circ_rna():
                circ_rna = traversal.circ_rna
                filtered_orfs = []
                for orf_i in cur_orfs:
                    orf_i_2 = orf_i.copy()
                    if not orf_i_2.is_valid_orf(out_node, self.subgraphs, circ_rna):
                        continue

                    if not orf_i_2.node_is_at_least_one_loop_downstream(
                                out_node, self.subgraphs, circ_rna):
                        for v in out_node.variants:
                            if v.not_cleavage_altering() \
                                    and orf_i_2.is_compatible_with_variant(v):
                                orf_i_2.start_gain.add(v.variant)
                    orf_i_2.cleavage_gain = set(cleavage_gain)
                    filtered_orfs.append(orf_i_2)
                if not filtered_orfs:
                    cur_in_cds = False
            else:
                filtered_orfs = cur_orfs

            cursor = PVGCursor(target_node, out_node, cur_in_cds, filtered_orfs,
                cur_cleavage_gain, finding_start_site)
            traversal.stage(target_node, out_node, cursor)

        for node, orfs, is_start_codon, additional_variants in node_list:
            if traversal.find_ass and any(o == traversal.known_orf_tx[0] for o in orfs):
                continue
            traversal.pool.add_miscleaved_sequences(
                node=node,
                orfs=orfs,
                cleavage_params=self.cleavage_params,
                check_variants=traversal.check_variants,
                is_start_codon=is_start_codon,
                additional_variants=additional_variants,
                denylist=self.denylist,
                leading_node=target_node,
                subgraphs=self.subgraphs,
                circ_rna=traversal.circ_rna,
                backsplicing_only=traversal.backsplicing_only
            )
        for node in trash:
            self.remove_node(node)

    def jsonfy(self):
        """ Create node and edge list from a PVG object. """
        queue = deque([self.root])
        node_index:dict[int, int] = {}
        cur_index = 0

        if self.has_known_orf():
            known_start_aa = int((self.known_orf[0] - self.known_reading_frame_index())/3)

        nodes = []
        while queue:
            cur = queue.pop()
            if id(cur) in node_index:
                continue
            if not cur.get_out_nodes():
                continue

            node = {
                'index': cur_index
            }
            node_index[id(cur)] = cur_index


            node['seq'] = str(cur.seq.seq) if cur.seq else ''
            variants = set()
            for v in cur.variants:
                if v.variant.attrs.get('MERGED_MNV'):
                    variants.update(v.variant.attrs['INDIVIDUAL_VARIANT_IDS'])
                else:
                    variants.add(v.variant.id)
            node['variants'] = list(variants)
            node['rf_index'] = cur.reading_frame_index

            if self.has_known_orf():
                if cur.seq:
                    i = cur.seq.get_query_index(
                        ref_index=known_start_aa,
                        seqname=self.id,
                        reading_frame=self.known_reading_frame_index()
                    )
                    node['has_start'] = i > -1
                else:
                    node['has_start'] = False
            else:
                node['has_start'] = cur.seq and 'M' in cur.seq.seq
            node['is_stop'] = cur.seq and cur.seq.seq == '*'
            nodes.append(node)

            for out_node in cur.get_out_nodes():
                if out_node.reading_frame_index == cur.reading_frame_index:
                    queue.append(out_node)
                else:
                    queue.appendleft(out_node)

            cur_index += 1

        edges = []
        queue = deque([self.root])
        visited:set[int] = set()
        while queue:
            cur = queue.pop()
            if id(cur) in visited:
                continue
            visited.add(id(cur))
            for out_node in cur.get_out_nodes():
                try:
                    source_node = node_index[id(cur)]
                    target_node = node_index[id(out_node)]
                except KeyError:
                    continue
                edge = {'source': source_node, 'target': target_node}
                edges.append(edge)
            for out_node in cur.get_out_nodes():
                queue.appendleft(out_node)

        return {
            'nodes': nodes,
            'edges': edges
        }


class PVGCursor():
    """ Helper class for cursors when graph traversal to call peptides. """
    def __init__(self, in_node:PVGNode, out_node:PVGNode, in_cds:bool,
            orfs:List[PVGOrf]=None, cleavage_gain:List[VariantRecord]=None,
            finding_start_site:bool=True):
        """ constructor """
        self.in_node = in_node
        self.out_node = out_node
        self.in_cds = in_cds
        self.cleavage_gain = cleavage_gain or []
        self.orfs = orfs or []
        self.finding_start_site = finding_start_site


class PVGTraversal():
    """ PVG Traversal. The purpose of this class is to facilitate the graph
    traversal to call variant peptides.
    """
    def __init__(self, check_variants:bool, check_orf:bool,
            pool:VariantPeptideDict, known_orf_tx:Tuple[int,int]=None,
            known_orf_aa:Tuple[int,int]=None, circ_rna:CircRNAModel=None,
            queue:Deque[PVGCursor]=None,
            stack:Dict[PVGNode, Dict[PVGNode, PVGCursor]]=None,
            orf_assignment:str='max', backsplicing_only:bool=False,
            find_ass:bool=False):
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

    def stage(self, in_node:PVGNode, out_node:PVGNode, cursor:PVGCursor):
        """ When a node is visited through a particular edge during the variant
        peptide finding graph traversal, it is staged until all inbond edges
        are visited. """
        in_nodes = self.stack.setdefault(out_node, {})

        if in_node in in_nodes:
            return
        in_nodes[in_node] = cursor

        if len(in_nodes) != len(out_node.in_nodes):
            return

        curs = list(self.stack[out_node].values())

        is_circ_rna = any(x.variant.is_circ_rna() for x in out_node.variants)

        if self.known_orf_aa[0] is not None:
            if out_node.reading_frame_index == self.known_reading_frame_index():
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
