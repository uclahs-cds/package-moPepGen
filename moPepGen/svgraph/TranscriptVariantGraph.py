""""""
from __future__ import annotations
from moPepGen.svgraph.PeptideVariantGraph import PeptideVariantGraph
from moPepGen.svgraph.VariantRecordWithCoordinate import VariantRecordWithCoordinate
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from typing import List, Set, Tuple, Dict
import copy
from collections import deque
from Bio import SeqUtils
from Bio.Seq import Seq, MutableSeq
from moPepGen import get_equivalent
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.dna.MatchedLocation import MatchedLocation
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord
from moPepGen.aa import AminoAcidSeqRecord
from moPepGen import svgraph


class TranscriptVariantGraph():
    """
    Defines the DAG(directed acyclic graph) data structure for the algorithm to
    find variant peptides raised by genomic variants.

    The algorithm has two main steps:

    1. Create the DAG with transcript variant records. For each variant, it
    first breaks the transcript into three nodes, as 'prior', 'between', and
    'post' to the variant location, and linked with edges that labeled as
    "reference". The variant is then linked to the pior and post node with
    edges labeled as 'variant_start' and 'variant_end'.

    2. Walk through the graph, find enzymatic sites, translate to peptides,
    and add to the final results if it contains a variant and does not have
    the same sequence as any canonical peptides from the reference transcript
    sequence.

    Attributes:
        nodes (Set[svgraph.Node]): Set of nodes. Each node ha
        edges (Set[svgraph.Edge]): Set of edges. Each 
        seq (DNASeqRecordWithCoordinates): The original sequence of the 
            transcript (reference).
        transcript_id (str)
        variants List[svgraph.VariantRecordWithCoordinate]: All variant
            records.
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates, transcript_id:str,
            protein:AminoAcidSeqRecord):
        """ Constructor to create a TranscriptVariantGraph object.
        
        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the 
                transcript (reference).
            transcript_id (str)
        """
        if seq.locations == []:
            seq.locations = [MatchedLocation(
                query=FeatureLocation(start=0, end=len(seq)),
                end=FeatureLocation(start=0, end=len(seq))
            )]
        self.seq = seq
        node = svgraph.Node(seq)
        self.root = node
        self.transcript_id = transcript_id
        self.variants = []
        self.protein = protein
    
    # def __repr__(self) -> str:
    #     """"""
    #     ref = MutableSeq('')
    #     cur = self.root
    #     while cur:
    #         if cur.seq is None:
    #             cur = cur.get_reference_next()
    #             continue
    #         if ref:
    #             ref.append('-')
    #         for i in cur.seq.seq:
    #             ref.append(i)
    #     out = deque([ref])
    #     up = True
    #     cur = self.root

    
    def variants_are_adjacent(self, edge1:svgraph.Edge, edge2:svgraph.Edge) -> bool:
        """ Deprecated If two edges are both variants. """
        if edge1 is None:
            return False
        if edge2 is None:
            return False
        return edge1.type == 'variant_start' and edge2.type == 'variant_start'
    
    def remove_edge(self, edge:svgraph.Edge) -> None:
        """ Removes an edge from the graph.
        
        The edge will be removed from
        the outbond edge list of its inbond node, and from the inbond edge
        list of its outbond node. The inbond or outbond node is removed if it
        became an orphan (do not have any edge connected) after the removal
        of the edge.
        
        Args:
            edge (svgraph.Edge): The edge to be removed.
        """
        if edge in edge.in_node.out_edges:
            edge.in_node.out_edges.remove(edge)
        if edge in edge.out_node.in_edges:
            edge.out_node.in_edges.remove(edge)
    
    def add_edge(self, in_node:svgraph.Node, out_node:svgraph.Node, type:str
            ) -> svgraph.Edge:
        """ Add an edge between two nodes.
        
        An edge is created abd added to the the outbond edge list of the inbond
        node, as well as the intbond edge list of the output node, then added
        to the edge list of the object.

        Args:
            in_node (svgraph.Node): the inbond node of the edge to add.
            out_node (svgraph.Node): the outbond node of the edge to add.
            type (str): the type of the edge to add. Can be either 'reference',
                'variant_start', and 'variant_end'.
        
        Returns:
            The new added edge.
        """
        edge = svgraph.Edge(in_node, out_node, type)
        in_node.out_edges.add(edge)
        out_node.in_edges.add(edge)
        return edge
    
    def remove_node(self, node:svgraph.Node):
        """"""
        while node.in_edges:
            edge = node.in_edges.pop()
            self.remove_edge(edge)

    def splice(self, node:svgraph.Node, i:int, type:str
            ) -> Tuple[svgraph.Node, svgraph.Node]:
        """ Splice a node into two separate nodes.

        The sequence of a given node is spliced at the given position. Two
        nodes with smaller sequences are then created with a edge pointing
        from the upstream to the downstream node. Any node connected two the 
        original node is transfered two the new nodes.

        Args:
            node (svgraph.Node): The node to be spliced.
            i (int): The position of sequence to splice.
            type (str): The type 
            type (str): the type of the edge to add bewen the two new nodes.
                Can be either 'reference', 'variant_start', and 'variant_end'.
        
        Returns:
            The two new nodes.
        """
        left_node = svgraph.Node(node.seq[:i])
        right_node = svgraph.Node(node.seq[i:])

        while node.in_edges:
            edge = node.in_edges.pop()
            self.add_edge(edge.in_node, left_node, type=edge.type)
            self.remove_edge(edge)
        
        while node.out_edges:
            edge = node.out_edges.pop()
            self.add_edge(right_node, edge.out_node, type=edge.type)
            self.remove_edge(edge)
        
        self.add_edge(left_node, right_node, type=type)
        
        return left_node, right_node
    
    
    def apply_variant(self, node:svgraph.Node, variant:VEPVariantRecord) -> svgraph.Node:
        """ Apply a given variant to the graph.
        
        For a given genomic variant with the coordinates of the transcript,
        the reference transcript sequence is first spliced into three fragments
        (nodes) as 'prior', 'between', and 'post' to the variant location.
        Unless the pior or post are already spliced. A new node with the
        sequence of the alt sequence is created and linked to the 'prior' and
        'post' nodes.

        Args:
            node [node]: The node where the variant to be add.
            variant [VEPVariantRecord]: The variant record.
        
        Returns:
            The reference node (unmutated) with the same location as the
            node is returned.
        """
        variant_start = variant.variant.location.start
        variant_end = variant.variant.location.end
        node_start = node.seq.locations[0].ref.start
        node_end = node.seq.locations[-1].ref.end
        if variant_start < node_start or variant_start > node_end:
            raise ValueError('Variant out of range')
        seq = DNASeqRecordWithCoordinates(
            seq=Seq(variant.variant.alt),
            locations=[],
            orf=None
        )
        variant_with_coordinates = svgraph.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=0, end=len(seq))
        )
        var_node = svgraph.Node(seq=seq, variants=[variant_with_coordinates])
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
            while cur.seq.locations[-1].ref.end <= variant_end:
                cur = cur.get_reference_next()
            index = cur.seq.get_query_index(variant_end)
            if index == 0:
                self.add_edge(var_node, cur, 'variant_end')
            else:
                left, right = self.splice(cur, index, 'reference')
                self.add_edge(var_node, right, 'variant_end')
        
        # returns the node with the first nucleotide of the input node.
        if head is not None:
            return head
        if mid is not None:
            return mid
        return body     
        
    
    def create_variant_graph(self, variants:List[VEPVariantRecord]) -> None:
        """ Create a variant graph.

        With a list of genomic variants, incorprate each variant into the
        graph.

        Args:
            variant [VEPVariantRecord]: The variant record.
        """
        self.variants += variants
        variant_iter = variants.__iter__()
        variant = next(variant_iter)
        cur = self.root

        while variant:
            
            if cur.seq.locations[0].ref.start > variant.variant.location.start:
                variant = next(variant_iter, None)
                continue
            # NOTE(CZ): a silent mutation may be not silent for cases that
            # a frameshifting mutation occured in the upstream!
            if variant.is_silent():
                variant = next(variant_iter, None)
                continue

            if cur.seq.locations[-1].ref.end <= variant.variant.location.start:
                try:
                    cur = cur.get_reference_next()
                except StopIteration:
                    break
                continue

            cur = self.apply_variant(cur, variant)
            if len(cur.in_edges) == 0:
                self.root = cur
            variant = next(variant_iter, None)
    

    def align_varints(self, node:svgraph.Node
            ) -> Tuple[svgraph.Node, List[svgraph.Node]]:
        r""" Aligns all variants at that overlaps to the same start and end
        position.

        For example:

                 T--                     TG--
                /   \                   /    \
            ATGG-TCT-G-CCCT   ->    ATGG-TCTG-CCCT
                    \ /                 \    /
                     A                   TCTA
        """
        start_node = node
        end_node = node.find_farthest_node_with_overlap()

        new_nodes = set()
        queue = deque()
        
        # start by removing the out edges from the start node.
        while start_node.out_edges:
            edge = start_node.out_edges.pop()
            out_node:svgraph.Node = edge.out_node
            self.remove_edge(edge)
            queue.appendleft(out_node)
            new_nodes.add(out_node)
        
        while queue:
            cur:svgraph.Node = queue.pop()
            if len(cur.out_edges) == 1 and \
                    next(iter(cur.out_edges)).out_node == end_node:
                continue
            new_nodes.remove(cur)
            
            while cur.out_edges:
                out_edge = cur.out_edges.pop()
                out_node:svgraph.Node = out_edge.out_node

                # If the next node has any variants, shift the location by the
                # length of the cur node's sequence.
                new_variants = []
                for variant in out_node.variants:
                    new_variants.append(variant.shift(len(cur.seq)))
            
                # create new node with the combined sequence
                new_node = svgraph.Node(
                    seq=cur.seq + out_node.seq,
                    variants=cur.variants + new_variants,
                    frameshifts=cur.framshifts + out_node.framshifts
                )
                for edge in out_node.out_edges:
                    edge_type = 'variant_end' if new_node.variants \
                        else 'reference'
                    self.add_edge(new_node, edge.out_node, type=edge_type)
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


    def expand_alignments(self, node:svgraph.Node
            ) -> Tuple[svgraph.Node, List[svgraph.Node]]:
        """ Expand the aligned variants into the range of codons. For
        frameshifting mutations, a copy of each downstream node will be
        created and branched out.

        Sequence after stop codons are trimmed off. However for frameshifting
        mutations where no stop codon is found, no trimming will be done. So
        the trailing sequence may not be multiple of 3. 
        
        For example:

             TG--                   GTG-CCCT   
            /    \                 /
        ATGG-TCTG-CCCT    ->    ATG-GTCTGC-CCT
            \    /                 \      /
             TCTA                   GTCTAC

        """
        branches = []
        trash_nodes = set()
        # how many nt we should give to the mutation.
        left_index = len(node.seq) - len(node.seq) % 3
        left_over = node.seq[left_index:]
        node.seq = node.seq[:left_index]
        ref_node = node.get_reference_next()
        original_ref_seq_len = len(ref_node.seq)
        end_node = ref_node.get_reference_next()
        # number of nt we should get from the reference node after mutation
        right_index = len(left_over) + len(ref_node.seq)
        right_index = (3 - right_index % 3) % 3
        if end_node:
            right_over = end_node.seq[:right_index]

        for edge in node.out_edges:
            out_node:svgraph.Node = edge.out_node

            out_node.seq = left_over + out_node.seq

            if out_node is ref_node:
                if end_node:
                    out_node.seq = out_node.seq + right_over
                continue

            variants = []
            for variant in out_node.variants:
                variants.append(variant.shift(len(left_over)))
            out_node.variants = variants

            # Look for stop codon. If stop codon is found, remove all
            # downstream
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

            # branch out framshifting mutations
            original_out_node_len = len(out_node.seq) - len(left_over)
            if (original_ref_seq_len - original_out_node_len) % 3 != 0:
                new_out_node:svgraph.Node = out_node.deepcopy()
                edge.out_node = new_out_node
                new_out_node.in_edges.add(edge)
                out_node.in_edges.remove(edge)
                
                # carry over nts to make codons if need
                new_next_node = new_out_node.get_reference_next()
                if len(new_out_node.seq) % 3 > 0 and new_next_node:
                    new_right_index = 3 - (len(new_out_node.seq) % 3)
                    new_out_node.seq = new_out_node.seq + \
                        new_next_node.seq[:new_right_index]
                    new_next_node.seq = new_next_node.seq[new_right_index:]

                more_frameshifts = [variant.variant for variant in \
                    new_out_node.variants]
                new_out_node.framshifts.extend(more_frameshifts)

                # the original frameshifting node should be then removed
                trash_nodes.add(out_node)

                new_cursor = svgraph.CursorForCodonFinding(
                    node=new_out_node,
                    position_to_orf='mid',
                    search_orf=True
                )
                branches.append(new_cursor)
                continue

            out_node.seq = out_node.seq + right_over
        
        for trash in trash_nodes:
            self.remove_node(trash)
        
        if end_node:
            end_node.seq = end_node.seq[right_index:]
            new_cursor = svgraph.CursorForCodonFinding(
                node=end_node,
                position_to_orf='mid'
            )
            return new_cursor, branches
        else:
            return None, branches
        
    def fit_into_codons(self) -> None:
        """ This takes the variant graph, and fit all cleave sites into the
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
        # Add a null nood as root. This might be moved to the
        # create_variant graph method.
        original_root = self.root
        new_root = svgraph.Node(None)
        self.root = new_root
        self.add_edge(self.root, original_root, 'reference')
        orf_start = self.seq.orf.start
        orf_end = self.seq.orf.end
        index = original_root.seq.get_query_index(orf_start)
        cur = svgraph.CursorForCodonFinding(
            node=original_root,
            position_to_orf='up',
            search_orf=False,
            start_codon_index=index
        )
        queue = deque([cur])

        while queue:
            cur:svgraph.CursorForCodonFinding = queue.pop()
            # not considering when there is mutation before the start codon
            # This should be the root
            if cur.position_to_orf == 'up' and \
                    cur.start_codon_index != -1:
                node:svgraph.Node = cur.node
                index = cur.start_codon_index
                node.seq = node.seq[index:]

                # When there is any start codon loss mutation, each mutation
                # will start it's own branch. The root will then become a null
                # node, and all branches will point to it.
                if len(node.seq) < 3:
                    edge:svgraph.Edge
                    for edge in cur.out_edges:
                        next_seq = cur.seq + edge.out_node.seq
                        variant = edge.out_node.variant
                        
                        # skip if this is a start retained mutation.
                        if variant and variant.is_nsv() and \
                                    next_seq.seq.startswith('ATG'):
                            self.remove_node(edge.out_node)
                            self.remove_edge(edge)
                            continue
                        
                        edge.out_node.seq = next_seq
                        edge.in_node = self.root

                        # reference, continue 
                        if variant is None and edge.type == 'reference':
                            new_cursor = svgraph.CursorForCodonFinding(
                                node=edge.out_node,
                                position_to_orf='mid',
                            )
                            queue.append(new_cursor)
                            continue
                        
                        # frameshifting mutation
                        edge.out_node = edge.out_node.deepcopy()
                        more_framshifts = set([variant.variant for variant in
                            edge.out_node.variants])
                        edge.out_node.framshifts.extend(more_framshifts)
                        new_cursor = svgraph.CursorForCodonFinding(
                            node=edge.out_node,
                            search_orf=True
                        )
                        queue.appendleft(edge.out_node)
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
                
                # no start codon found in the current cursor
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
                    start_codon=index
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
                    node = self.align_varints(node)
                main, branches = self.expand_alignments(node)
                for branch in branches:
                    queue.appendleft(branch)
                if main:
                    queue.append(main)
                continue
        return
    
    def translate(self, transcript_id:str, protein_id:str, gene_id:str
            ) -> svgraph.PeptideVariantGraph:
        """ Converts a DNA transcript variant graph into a peptide variant
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
                out_node:svgraph.Node = edge.out_node
                if out_node in visited:
                    out_pnode = visited[out_node]
                    pnode.add_out_edge(out_pnode)
                    continue

                seq = out_node.seq.translate()
                stop_index = seq.seq.find('*')
                
                if stop_index > -1:
                    seq = seq[:stop_index]
                
                seq.transcript_id = transcript_id
                seq.protein_id = protein_id
                seq.gene_id = gene_id

                # translate the dna variant location to peptide coordinates.
                variants = []
                for variant in out_node.variants:
                    start = int(variant.location.start/3)
                    end = int(variant.location.end/3)
                    if start > len(seq):
                        continue
                    new_variant = VariantRecordWithCoordinate(
                        variant=variant.variant,
                        location=FeatureLocation(start=start, end=end)
                    )
                    variants.append(new_variant)
                    
                new_pnode = svgraph.PeptideNode(
                    seq=seq,
                    variants=variants,
                    frameshifts=set(out_node.framshifts)
                )
                pnode.add_out_edge(new_pnode)
                visited[out_node] = new_pnode
                if stop_index > -1:
                    pgraph.add_stop(new_pnode)
                else:
                    queue.appendleft((out_node, new_pnode))
        return pgraph
    
    def walk_and_splice(self, rule:str, exception:str=None, miscleavage:int=2,
            min_mw:float=500.) -> Set[AminoAcidSeqRecord]:
        """ Walk through the graph and finds variant peptides.

        The part of the sequence before start codon is first skipped. It then
        walks through each path of the graph. For each iteration, it only looks
        at the sequence with at least three enzymatic cleave sites. Piptides
        with specified number of miscleavages are extracted, and added to the
        final report set if not seen in the carnonincal peptide pool of the 
        transcript's reference sequence.

        For the case of a frameshifting mutation, the start codon and/or stop
        codon is searched while walking through the path.

        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            miscleavage (int): Number of miscleavage to allow. Defaults to 2.
            min_mw (float): Minimal molecular weight of the peptides to report.
                Defaults to 500.

        Returns:
            A set of variant peptides.
        """
        ref_orf = self.root.seq.orf

        protein = self.protein
        if protein.seq.startswith('X'):
            protein.seq = protein.seq.lstrip('X')

        carnonical_peptides = protein.enzymatic_cleave(
            rule=rule,
            miscleavage=miscleavage,
            exception=exception,
            min_mw=min_mw
        )
        carnonical_peptides = set([peptide.seq for peptide in \
            carnonical_peptides])
        # DNASeqRecordWithCoordinates only hashes the sequence, so here we are
        # only keeping the unique sequences of the carnonical peptides,
        # regardless of the location.
        carnonical_peptides = set(carnonical_peptides)
        incarnonical_peptides = set()

        queue = deque()
        location = FeatureLocation(start=0, end=len(self.root.seq))
        node_with_location = NodeLocation(node=self.root, location=location)
        queue.appendleft(Cursor(NodePath([node_with_location])))
        visited_sequences = set()
        visited_nodes = set()

        counter_for_debuging = 0
        while queue:
            counter_for_debuging += 1
            if counter_for_debuging % 500000 == 0:
                counter_for_debuging
                pass
            cursor:Cursor = queue.pop()

            if cursor.position_to_orf == 'up':
                # when it passes the ref_orf start site, and found that there
                # is a start codon lost mutation.
                if cursor.search_orf:
                    orf_start = cursor.seq.find_orf_first()
                    for edge in cursor.nodes[-1].node.out_edges:
                        node = edge.out_node
                        if orf_start == -1:
                            # It is truncated because the first part is already
                            # check that there isn't a start codon.
                            new_cursor = cursor.join(node=node, start=-8,
                                search_orf=True)
                        else:
                            new_cursor = cursor.join(node=node,
                                start=orf_start, position_to_orf='mid',
                                search_orf=True)
                        queue.appendleft(new_cursor)
                        continue
                # when it is before the ref_orf start
                if cursor.seq.locations[-1].ref.end - 3 < ref_orf.start:
                    for edge in cursor.nodes[-1].node.out_edges:
                        new_cursor = cursor.join(node=edge.out_node, start=-2)
                        queue.appendleft(new_cursor)
                    continue
                # when it is after the ref_orf start, we then need to check if
                # start codon is indeed there
                for edge in cursor.nodes[-1].node.out_edges:
                    index = cursor.seq.get_query_index(ref_orf.start)
                    new_cursor = cursor.join(edge.out_node, index, 'mid')
                    if new_cursor.seq.seq[:3] != 'ATG':
                        new_cursor.search_orf = True
                        new_cursor.position_to_orf = 'up'
                    queue.appendleft(new_cursor)
                continue
            
            # If a frameshifting mutation occured after start codon, the
            # stop codon needs to be searched.
            if cursor.frameshifted and not cursor.search_orf \
                    and cursor.position_to_orf != 'down':
                cursor.search_orf = True
            
            # In case of frameshifting mutation, that the ORF needs to be
            # searched
            if cursor.search_orf:
                stop_codon = cursor.seq.find_stop_codon()
                if stop_codon > -1 :
                    cursor = cursor[:stop_codon]
                    cursor.search_orf = False
                    cursor.position_to_orf = 'down'
                if cursor.nodes[-1].node.out_edges is None:
                    cursor.position_to_orf = 'down'
            
            if cursor.position_to_orf == 'mid' and \
                    len(cursor.seq.locations) > 0 and \
                    cursor.seq.locations[-1].ref.end > ref_orf.end:
                # means the stop codon is in the current cursor
                stop_codon_index = cursor.seq\
                    .get_query_index(ref_orf.end)
                new_cursor = cursor[:stop_codon_index]
                new_cursor.position_to_orf = 'down'
                queue.appendleft(new_cursor)
                continue
            
            if len(cursor.cleave_sites) == 0:
                cursor.cleave_sites = cursor.seq\
                    .find_all_enzymatic_cleave_sites(
                        rule=rule, exception=exception)
            
            # if the seq is too short
            if len(cursor.cleave_sites) < miscleavage + 1 and \
                    cursor.position_to_orf == 'mid':
                for edge in cursor.nodes[-1].node.out_edges:
                    new_cursor = cursor.join(edge.out_node)
                    new_cursor.search_more_cleave_sites(rule=rule,
                            exception=exception)
                    queue.appendleft(new_cursor)
                continue
            
            if not cursor.has_any_variants():
                for edge in cursor.nodes[-1].node.out_edges:
                    if len(cursor.cleave_sites) > miscleavage:
                        start = cursor.cleave_sites[-(miscleavage + 1)]
                    else:
                        start = 0
                    new_cursor = cursor.join(edge.out_node, start)
                    queue.appendleft(new_cursor)
                continue
            
            visited_len_before = len(visited_nodes)
            visited_nodes.add(cursor.nodes)
            visited_len_after = len(visited_nodes)
            if visited_len_before == visited_len_after:
                continue

            cleave_iter = iter(cursor.cleave_sites)

            hit_stop_codon = False
            try:
                cleave_site = next(cleave_iter)
            except StopIteration:
                cleave_site = len(cursor.seq)
                hit_stop_codon = True
            first_cleave_site = cleave_site
            i = 0
            while i <= miscleavage:
                i += 1
                if cleave_site is None:
                    hit_stop_codon = True
                    break
                seq = cursor.seq[:cleave_site]
                
                # check if the sequence is already visited
                visited_len_before = len(visited_sequences)
                seq_data_to_visit = str(seq.seq)
                visited_sequences.add(seq_data_to_visit)
                visited_len_after = len(visited_sequences)
                if visited_len_after == visited_len_before:
                    cleave_site = next(cleave_iter, None)
                    continue
                
                # remove peptides without variants
                variants = cursor.get_variants_between(0, cleave_site)
                if len(variants) == 0:
                    cleave_site = next(cleave_iter, None)
                    continue
                
                peptide = seq.translate(to_stop=True)
                # means there may be a stop codon gained mutation
                # TODO(CZ): valiate the stop is caused by mutation
                # remove peptides with small MW
                mw = SeqUtils.molecular_weight(peptide.seq, 'protein')
                if mw < min_mw:
                    cleave_site = next(cleave_iter, None)
                    continue

                # check if peptide is seen in WT peptide
                if peptide.seq in carnonical_peptides:
                    cleave_site = next(cleave_iter, None)
                    continue
                
                # If the variant peptide is already reported, add the variant
                # information to that record, instead of reporting it again.
                variant_info = ''
                for variant in variants:
                    variant_info += ('|' + str(variant))
                same_peptide = get_equivalent(incarnonical_peptides,
                    peptide)
                if same_peptide:
                    same_peptide.id += variant_info
                    same_peptide.name = same_peptide.id
                    same_peptide.description = same_peptide.id
                    cleave_site = next(cleave_iter, None)
                    continue
                peptide.id = self.transcript_id + variant_info
                peptide.name = peptide.id
                peptide.description = peptide.id
                incarnonical_peptides.add(peptide)
            
            # If we are hitting the stop codon, we are done with this path.
            if hit_stop_codon:
                continue
            new_cursor = cursor[first_cleave_site:]
            queue.appendleft(new_cursor)
        return incarnonical_peptides