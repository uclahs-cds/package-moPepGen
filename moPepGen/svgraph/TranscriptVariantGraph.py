""" Module for transcript (DNA) variant graph """
from __future__ import annotations
from typing import List, Tuple, Deque, Set
from collections import deque
import copy
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import dna, seqvar, svgraph


class TranscriptVariantGraph():
    """ Defines the DAG data structure for the transcript and its variants.

    Attributes:
        nodes (Set[svgraph.DNANode]): Set of nodes. Each node ha
        edges (Set[svgraph.DNAEdge]): Set of edges. Each
        seq (DNASeqRecordWithCoordinates): The original sequence of the
            transcript (reference).
        transcript_id (str)
        variants List[seqvar.VariantRecordWithCoordinate]: All variant
            records.
    """
    def __init__(self, seq:dna.DNASeqRecordWithCoordinates|None,
            transcript_id:str):
        """ Constructor to create a TranscriptVariantGraph object.

        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the
                transcript (reference).
            transcript_id (str)
        """
        if seq:
            if seq.locations == []:
                seq.locations = [dna.MatchedLocation(
                    query=FeatureLocation(start=0, end=len(seq)),
                    ref=FeatureLocation(start=0, end=len(seq))
                )]
        self.seq = seq
        node = svgraph.DNANode(seq)
        self.root = node
        self.transcript_id = transcript_id
        self.variants = []

    @staticmethod
    def variants_are_adjacent(edge1:svgraph.DNAEdge, edge2:svgraph.DNAEdge
            ) -> bool:
        """ Deprecated If two edges are both variants. """
        if edge1 is None:
            return False
        if edge2 is None:
            return False
        return edge1.type == 'variant_start' and edge2.type == 'variant_start'

    @staticmethod
    def remove_edge(edge:svgraph.DNAEdge) -> None:
        """ Removes an edge from the graph.

        The edge will be removed from
        the outbond edge list of its inbond node, and from the inbond edge
        list of its outbond node. The inbond or outbond node is removed if it
        became an orphan (do not have any edge connected) after the removal
        of the edge.

        Args:
            edge (svgraph.DNAEdge): The edge to be removed.
        """
        if edge in edge.in_node.out_edges:
            edge.in_node.out_edges.remove(edge)
        if edge in edge.out_node.in_edges:
            edge.out_node.in_edges.remove(edge)

    @staticmethod
    def add_edge(in_node:svgraph.DNANode, out_node:svgraph.DNANode,
            _type:str) -> svgraph.DNAEdge:
        """ Add an edge between two nodes.

        An edge is created abd added to the the outbond edge list of the inbond
        node, as well as the intbond edge list of the output node, then added
        to the edge list of the object.

        Args:
            in_node (svgraph.DNANode): the inbond node of the edge to add.
            out_node (svgraph.DNANode): the outbond node of the edge to add.
            type (str): the type of the edge to add. Can be either 'reference',
                'variant_start', and 'variant_end'.

        Returns:
            The new added edge.
        """
        edge = svgraph.DNAEdge(in_node, out_node, _type)
        in_node.out_edges.add(edge)
        out_node.in_edges.add(edge)
        return edge

    def remove_node(self, node:svgraph.DNANode):
        """ Remove a node from its inboud and outbound node. """
        while node.in_edges:
            edge = node.in_edges.pop()
            self.remove_edge(edge)
        while node.out_edges:
            edge = node.out_edges.pop()
            self.remove_edge(edge)

    def add_null_root(self):
        """ Adds a null node to the root. """
        original_root = self.root
        new_root = svgraph.DNANode(None)
        self.root = new_root
        self.add_edge(self.root, original_root, 'reference')

    def splice(self, node:svgraph.DNANode, i:int, _type:str
            ) -> Tuple[svgraph.DNANode, svgraph.DNANode]:
        """ Splice a node into two separate nodes.

        The sequence of a given node is spliced at the given position. Two
        nodes with smaller sequences are then created with a edge pointing
        from the upstream to the downstream node. Any node connected two the
        original node is transfered two the new nodes.

        Args:
            node (svgraph.DNANode): The node to be spliced.
            i (int): The position of sequence to splice.
            type (str): The type
            type (str): the type of the edge to add bewen the two new nodes.
                Can be either 'reference', 'variant_start', and 'variant_end'.

        Returns:
            The two new nodes.
        """
        left_node = svgraph.DNANode(node.seq[:i])
        right_node = svgraph.DNANode(node.seq[i:])

        while node.in_edges:
            edge = node.in_edges.pop()
            self.add_edge(edge.in_node, left_node, _type=edge.type)
            self.remove_edge(edge)

        while node.out_edges:
            edge = node.out_edges.pop()
            self.add_edge(right_node, edge.out_node, _type=edge.type)
            self.remove_edge(edge)

        self.add_edge(left_node, right_node, _type=_type)

        if self.root is node:
            self.root = left_node

        return left_node, right_node


    def apply_variant(self, node:svgraph.DNANode, variant:seqvar.VariantRecord
            ) -> svgraph.DNANode:
        """ Apply a given variant to the graph.

        For a given genomic variant with the coordinates of the transcript,
        the reference transcript sequence is first spliced into three fragments
        (nodes) as 'prior', 'between', and 'post' to the variant location.
        Unless the pior or post are already spliced. A new node with the
        sequence of the alt sequence is created and linked to the 'prior' and
        'post' nodes.

        Args:
            node [node]: The node where the variant to be add.
            variant [seqvar.VariantRecord]: The variant record.

        Returns:
            The reference node (unmutated) with the same location as the
            node is returned.
        """
        variant_start = variant.location.start
        variant_end = variant.location.end
        node_start = node.seq.locations[0].ref.start
        node_end = node.seq.locations[-1].ref.end
        if variant_start < node_start or variant_start > node_end:
            raise ValueError('Variant out of range')
        seq = dna.DNASeqRecordWithCoordinates(
            seq=Seq(variant.alt),
            locations=[],
            orf=None
        )
        variant_with_coordinates = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=0, end=len(seq))
        )
        var_node = svgraph.DNANode(
            seq=seq,
            variants=[variant_with_coordinates]
        )
        head = None
        mid = None
        # variant start
        if variant_start == node_start:
            # No need to splice, just add a edge to the prev node
            # This is the case that another variant happend at the same
            # location
            prev = node.get_reference_prev()
            self.add_edge(prev, var_node, 'variant_start')
            body = node
        else:
            index = node.seq.get_query_index(variant_start)
            head, body = self.splice(node, index, 'reference')
            self.add_edge(head, var_node, 'variant_start')

        # varaint end
        if variant_end < node_end:
            index = body.seq.get_query_index(variant_end)
            mid, tail = self.splice(body, index, 'reference')
            self.add_edge(var_node, tail, 'variant_end')
        else:
            # This is the case that the range of this variant is larger than
            # the node, such as deletion.
            cur = body
            while cur.seq.locations[-1].ref.end <= variant_end and \
                    cur.out_edges:
                cur = cur.get_reference_next()

            # this is to avoid the case that the mutation is at the very last
            # nt
            if cur.seq.locations[-1].ref.end > variant_end:
                index = cur.seq.get_query_index(variant_end)
                if index == 0:
                    self.add_edge(var_node, cur, 'variant_end')
                else:
                    _, right = self.splice(cur, index, 'reference')
                    self.add_edge(var_node, right, 'variant_end')

        # returns the node with the first nucleotide of the input node.
        if head is not None:
            return head
        if mid is not None:
            return mid
        return body


    def create_variant_graph(self, variants:List[seqvar.VariantRecord]
            ) -> None:
        """ Create a variant graph.

        With a list of genomic variants, incorprate each variant into the
        graph.

        Args:
            variant [seqvar.VariantRecord]: The variant record.
        """
        # Add a null nood as root.
        self.add_null_root()
        self.variants += variants
        variant_iter = variants.__iter__()
        variant = next(variant_iter)
        cur = self.root.get_reference_next()

        while variant and cur:

            if cur.seq.locations[0].ref.start > variant.location.start:
                variant = next(variant_iter, None)
                continue

            if cur.seq.locations[-1].ref.end <= variant.location.start:
                cur = cur.get_reference_next()
                continue

            cur = self.apply_variant(cur, variant)
            if len(cur.in_edges) == 0:
                self.root = cur
            variant = next(variant_iter, None)


    def create_branch(self, node:svgraph.DNANode) -> svgraph.DNANode:
        """ Create a branch by making a deep copy of the specified node and
        all its children. The branch attribute of the copied node is set
        to True. """
        new_node:svgraph.DNANode = node.deepcopy()
        edges = copy.copy(node.in_edges)
        while edges:
            edge = edges.pop()
            self.add_edge(edge.in_node, new_node, edge.type)
            self.remove_edge(edge)
        new_node.branch = True
        return new_node

    def merge_with_outbonds(self, node:svgraph.DNANode) -> List[svgraph.DNANode]:
        """ For a given node, merge it with all its outbound nodes. """

        in_edges = copy.copy(node.in_edges)
        out_edges:Set[svgraph.DNAEdge] = copy.copy(node.out_edges)
        out_nodes = []
        while out_edges:
            edge = out_edges.pop()
            out_node = edge.out_node
            out_node.seq = node.seq + out_node.seq

            for i,_ in enumerate(out_node.variants):
                out_node.variants[i] = out_node.variants[i].shift(len(node.seq))

            for in_edge in in_edges:
                self.add_edge(in_edge.in_node, out_node, edge.type)

            out_nodes.append(out_node)

        self.remove_node(node)

        return out_nodes

    def align_variants(self, node:svgraph.DNANode, branch_out_size:int=100
            ) -> svgraph.DNANode:
        r""" Aligns all variants at that overlaps to the same start and end
        position. Frameshifting mutations will be brached out

        For example:

                 T--                     T-G-CCCT
                /   \                   /
            ATGG-TCT-G-CCCT   ->    ATGG-TCTG-CCCT
                    \ /                 \    /
                     A                   TCTA

        Args:
            node (svgraph.DNANode): The node of which the outbound nodes will
                be aligned.
            branch_out_size (int): The size limit of the variant to forch
                creating a branch in th graph.
        Returns:
            The original input node.
        """
        start_node = node
        end_node = node.find_farthest_node_with_overlap()
        node_to_branch = start_node.next_node_to_branch_out(
            to_node=end_node,
            branch_out_size=branch_out_size
        )
        while node_to_branch:
            self.create_branch(node_to_branch)
            end_node = node.find_farthest_node_with_overlap()
            node_to_branch = start_node.next_node_to_branch_out(
                to_node=end_node,
                branch_out_size=branch_out_size
            )

        new_nodes = set()
        queue = deque()

        # start by removing the out edges from the start node.
        while start_node.out_edges:
            edge = start_node.out_edges.pop()
            out_node:svgraph.DNANode = edge.out_node
            self.remove_edge(edge)
            queue.appendleft(out_node)
            new_nodes.add(out_node)

        while queue:
            cur:svgraph.DNANode = queue.pop()
            if cur is end_node or ((len(cur.out_edges) == 1 and \
                    next(iter(cur.out_edges)).out_node == end_node)) or \
                    (not cur.out_edges):
                continue
            if cur.branch:
                continue

            # because the new node will be reconstructed
            new_nodes.remove(cur)

            while cur.out_edges:
                out_edge = cur.out_edges.pop()
                out_node:svgraph.DNANode = out_edge.out_node

                # If the next node has any variants, shift the location by the
                # length of the cur node's sequence.
                new_variants = []
                for variant in out_node.variants:
                    new_variants.append(variant.shift(len(cur.seq)))

                # create new node with the combined sequence
                frameshifts = cur.frameshifts
                frameshifts.update(out_node.frameshifts)
                new_node = svgraph.DNANode(
                    seq=cur.seq + out_node.seq,
                    variants=cur.variants + new_variants,
                    frameshifts=frameshifts,
                    branch=out_node.branch
                )
                for edge in out_node.out_edges:
                    edge_type = 'variant_end' if new_node.variants \
                        else 'reference'
                    self.add_edge(new_node, edge.out_node, _type=edge_type)
                queue.appendleft(new_node)
                new_nodes.add(new_node)
            # now remove the cur node from graph
            self.remove_node(cur)

        # add now nodes to the graph
        for new_node in new_nodes:
            if new_node.variants:
                in_edge_type = 'variant_start'
            else:
                in_edge_type = 'reference'
            self.add_edge(start_node, new_node, in_edge_type)

        return start_node

    def expand_alignments(self, node:svgraph.DNANode
            ) -> Tuple[svgraph.CursorForCodonFinding,
                    List[svgraph.CursorForCodonFinding]]:
        r""" Expand the aligned variants into the range of codons. For
        frameshifting mutations, a copy of each downstream node will be
        created and branched out.

        Sequence after stop codons are trimmed off. However for frameshifting
        mutations where no stop codon is found, no trimming will be done. So
        the trailing sequence may not be multiple of 3.

        For example:

             T-GCCCT                GTG-CCCT
            /                      /
        ATGG-TCTG-CCCT    ->    ATG-GTCTGC-CCT
            \    /                 \      /
             TCTA                   GTCTAC
        """
        branches = []
        # how many nt we should give to the mutation.
        if node.seq:
            left_index = len(node.seq) - len(node.seq) % 3
            left_over = node.seq[left_index:]
            node.seq = node.seq[:left_index]
        else:
            left_over = dna.DNASeqRecordWithCoordinates(Seq(''), [])

        ref_node = node.get_reference_next()
        # original_ref_seq_len = len(ref_node.seq)
        end_node = ref_node.get_reference_next()
        # number of nt we should get from the reference node after mutation
        right_index = len(left_over) + len(ref_node.seq)
        right_index = (3 - right_index % 3) % 3
        if end_node:
            right_over = end_node.seq[:right_index]

        needs_expand_forward = False
        if len(ref_node.seq.seq) + len(left_over) < 3:
            needs_expand_forward = True
        for edge in node.out_edges:
            for variant in edge.out_node.variants:
                if not variant.variant.is_frameshifting():
                    needs_expand_forward = True
                    break

        for edge in node.out_edges:
            out_node:svgraph.DNANode = edge.out_node

            out_node.seq = left_over + out_node.seq

            if out_node is ref_node:
                if end_node and needs_expand_forward:
                    out_node.seq = out_node.seq + right_over
                continue

            variants = []
            for variant in out_node.variants:
                variants.append(variant.shift(len(left_over)))
            out_node.variants = variants

            # Look for stop codon. If stop codon is found, remove all
            # downstream edges/nodes
            stop_codon_index = out_node.seq.find_stop_codon()
            if stop_codon_index > -1:
                out_node.seq = out_node.seq[:stop_codon_index + 3]
                # only keep the variants before the end of staop codon
                variants = [variant for variant in variants if
                    variant.location.start < stop_codon_index + 3]
                out_node.variants = variants
                while out_node.out_edges:
                    out_edge = out_node.out_edges.pop()
                    self.remove_edge(out_edge)
                continue

            if out_node.branch:
                # carry over nts to make codons if need
                branch_next_node = out_node.get_reference_next()
                if len(out_node.seq) % 3 > 0 and branch_next_node:
                    branch_right_index = 3 - (len(out_node.seq) % 3)
                    #
                    if len(branch_next_node.seq.seq) <= branch_right_index:
                        downstreams = self.merge_with_outbonds(branch_next_node)
                        # checks if any of the outbound nodes is still too short
                        out_node_needs_merge_with_its_outnodes = False
                        if len(downstreams) > 1 or \
                                len(downstreams[0].seq) <= branch_right_index:
                            out_node_needs_merge_with_its_outnodes = True

                        if out_node_needs_merge_with_its_outnodes:
                            new_out_nodes = self.merge_with_outbonds(out_node)
                            branch_next_node = new_out_nodes[0]\
                                .get_reference_next()
                            branches.append(branch_next_node)
                            continue
                        branch_next_node = downstreams[0]
                    out_node.seq = out_node.seq + \
                        branch_next_node.seq[:branch_right_index]
                    branch_next_node.seq = branch_next_node.seq[branch_right_index:]

                branches.append(branch_next_node)
                continue

            if end_node:
                out_node.seq = out_node.seq + right_over

        if needs_expand_forward:
            if end_node:
                return end_node, branches
            return None, branches
        return ref_node, branches

    def fit_into_codons(self, branch_out_size:int=100) -> None:
        r""" This takes the variant graph, and fit all cleave sites into the
        range of codons. For frameshifting mutations, the downstream part
        will be branched out.

        The case where there are mutations before the start codon isn't
        considered

        Example:

                 T--                    GTG-CCCT
                /   \                  /
            ATGG-TCT-G-CCCT   ->    ATG-GTCTGC-CCT
                    \ /                \      /
                     A                  GTCTAC
        """
        cur = svgraph.CursorForCodonFinding(
            node=self.root,
            position_to_orf='up',
            search_orf=False
        )
        queue = deque([cur])

        while queue:
            cur:svgraph.CursorForCodonFinding = queue.pop()
            if cur.node is self.root:
                if len(cur.node.out_edges) == 1:
                    main = cur.node.get_reference_next()
                else:
                    node = self.align_variants(cur.node, branch_out_size)
                    main,branches = self.expand_alignments(node)
                    queue.appendleft(main)
                    for branch in branches:
                        new_cursor = svgraph.CursorForCodonFinding(
                            node=branch,
                            position_to_orf='up',
                            search_orf=True
                        )
                        queue.appendleft(new_cursor)

                if self.seq.orf:
                    start_codon_index = main.seq\
                        .get_query_index(self.seq.orf.start)
                    search_orf = False
                else:
                    start_codon_index = None
                    search_orf = True
                new_cursor = svgraph.CursorForCodonFinding(
                    node=main,
                    position_to_orf='up',
                    search_orf=search_orf,
                    start_codon_index=start_codon_index
                )
                queue.append(new_cursor)
                continue

            # not considering when there is mutation before the start codon
            # This should be the root
            if cur.position_to_orf == 'up' and cur.start_codon_index != -1:
                node:svgraph.DNANode = cur.node
                index = cur.start_codon_index
                node.seq = node.seq[index:]

                # When there is any start codon loss mutation, each mutation
                # will start it's own branch.
                if len(node.seq) < 3:
                    edges = copy.copy(cur.node.out_edges)
                    while edges:
                        edge:svgraph.DNAEdge = edges.pop()
                        next_seq = cur.node.seq + edge.out_node.seq

                        # reference, continue
                        if not edge.out_node.variants:
                            new_cursor = svgraph.CursorForCodonFinding(
                                node=edge.out_node,
                                position_to_orf='mid',
                            )
                            queue.append(new_cursor)
                            continue

                        variant = edge.out_node.variants[0].variant

                        # skip if this is a start retained mutation.
                        if not variant.is_deletion() and \
                                next_seq.seq.endswith('ATG'):
                            self.remove_node(edge.out_node)
                            self.remove_edge(edge)
                            continue

                        edge.out_node.seq = next_seq
                        edge.in_node = self.root

                        # frameshifting mutation
                        out_node_copy = edge.out_node.deepcopy()
                        self.add_edge(edge.in_node, out_node_copy, edge.type)
                        self.remove_node(edge.out_node)
                        self.remove_edge(edge)
                        new_cursor = svgraph.CursorForCodonFinding(
                            node=out_node_copy,
                            search_orf=True
                        )
                        queue.appendleft(new_cursor)
                    continue

                # otherwise, there isn't any mutation at the start codon.
                # Continue to the next iteration
                cur.position_to_orf = 'mid'
                queue.append(cur)
                continue

            if cur.search_orf and cur.position_to_orf == 'up':
                # look for start codon
                node = cur.node
                index = node.seq.seq.find('ATG')

                # for start lost mutation or non-translatable sequences, if
                # start codon is not found in the current cursor, merge the
                # last nucleotides with each of the out nodes, and research.
                # This may be improved by skipping visited nodes.
                if index == -1:
                    for edge in node.out_edges:
                        out_node = edge.node
                        out_node.seq = node.seq[-2:] + out_node.seq
                        edge.in_node = node.in_node
                        self.remove_node(node)
                        new_cursor = svgraph.CursorForCodonFinding(
                            node=out_node,
                            search_orf=True
                        )
                        queue.appendleft(new_cursor)
                    continue

                # found start codon
                new_cursor = svgraph.CursorForCodonFinding(
                    node=node,
                    position_to_orf='up',
                    search_orf=True,
                    start_codon_index=index
                )
                queue.append(new_cursor)

            if cur.position_to_orf == 'mid':
                node = cur.node
                stop_codon_index = node.seq.find_stop_codon()
                if stop_codon_index != -1:
                    node.seq = node.seq[:stop_codon_index + 3]
                    while node.out_edges:
                        edge = node.out_edges.pop()
                        self.remove_edge(edge)
                    continue
                if not node.out_edges:
                    continue
                if len(node.out_edges) > 1:
                    node = self.align_variants(node, branch_out_size)
                main, branches = self.expand_alignments(node)

                if main:
                    new_cursor = svgraph.CursorForCodonFinding(
                        node=main,
                        position_to_orf='mid',
                        search_orf=False
                    )
                    queue.appendleft(new_cursor)
                for branch in branches:
                    new_cursor = svgraph.CursorForCodonFinding(
                        node=branch,
                        position_to_orf='mid',
                        search_orf=False
                    )
                    queue.appendleft(new_cursor)
                continue

    def translate(self) -> svgraph.PeptideVariantGraph:
        r""" Converts a DNA transcript variant graph into a peptide variant
        graph. A stop * is added to the end of all branches.

                GTG-CCCT              V-P-*
               /                     /
            ATG-GTCTGC-CCT    ->    M-VC-P-*
               \      /              \  /
                GTCTAC                VY
        """
        root = svgraph.PeptideNode(None)
        pgraph = svgraph.PeptideVariantGraph(root)

        queue = deque([(self.root, root)])
        visited = dict()

        while queue:
            dnode, pnode = queue.pop()
            if not dnode.out_edges:
                pgraph.add_stop(pnode)
                continue
            for edge in dnode.out_edges:
                out_node:svgraph.DNANode = edge.out_node
                if out_node in visited:
                    out_pnode = visited[out_node]
                    pnode.add_out_edge(out_pnode)
                    continue

                seq = out_node.seq.translate()
                stop_index = seq.seq.find('*')

                if stop_index > -1:
                    seq = seq[:stop_index]

                seq.transcript_id = self.transcript_id

                # translate the dna variant location to peptide coordinates.
                variants = []
                for variant in out_node.variants:
                    start = int(variant.location.start/3)
                    end = int(variant.location.end/3)
                    if start > len(seq):
                        continue
                    new_variant = seqvar.VariantRecordWithCoordinate(
                        variant=variant.variant,
                        location=FeatureLocation(start=start, end=end)
                    )
                    variants.append(new_variant)

                new_pnode = svgraph.PeptideNode(
                    seq=seq,
                    variants=variants,
                    frameshifts=set(out_node.frameshifts)
                )
                pnode.add_out_edge(new_pnode)
                visited[out_node] = new_pnode
                if stop_index > -1:
                    pgraph.add_stop(new_pnode)
                else:
                    queue.appendleft((out_node, new_pnode))
        return pgraph
