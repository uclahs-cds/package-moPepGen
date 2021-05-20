""" The DAG data structure and algorithm for finding variant peptides """
from __future__ import annotations
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from typing import List, Set, Tuple
from collections import deque
from Bio import SeqUtils
from moPepGen import get_equivalent
from moPepGen.dna import DNASeqRecordWithCoordinates
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.dna.MatchedLocation import MatchedLocation
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord


class Edge():
    """ Defines the edges in the TranscriptVariantGraph

    Attributes:
        in_node (Node): The inbond node.
        out_node (Node): The outbond node.
        type (str): The edge type. Must be either of orf_start, orf_end,
            variant_start, variant_end, cleave, or reference
    """
    def __init__(self, in_node:Node, out_node:Node,
            type:str):
        """ Constructor for Edge

        Args:
            in_node (Node): The inbond node.
            out_node (Node): The outbond node.
            type (str): The edge type. Must be either of orf_start, orf_end,
                mutation_start, mutation_end, cleave, or reference
        """
        self.in_node = in_node
        self.out_node = out_node
        edge_types = ['variant_start', 'variant_end', 'reference']
        if type not in edge_types:
            raise ValueError(f'type {type} not from {edge_types}')
        self.type = type

class Node():
    """ Defines the nodes in the TranscriptVariantGraph
    
    Attributes:
        in_edges (Set[Edge]): The inbonding edges.
        out_edges (Set[Edge]): The outbonding edges.
        seq (DNASeqRecord): The sequence.
        variant (VEPVariantRecord | None): The variant record or None for
            reference.
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates,
            variant:VEPVariantRecord=None):
        """ Constructor for Node.
        
        Args:
            seq (DNASeqRecord): The sequence.
            variant (VEPVariantRecord | None): The variant record or None for
                reference.
        """
        self.seq = seq
        self.in_edges = set()
        self.out_edges = set()
        self.variant = variant
    
    def __hash__(self):
        """ hash """
        locations = [(it.ref.start, it.ref.end) for it in self.seq.locations]
        return hash((str(self.seq.seq), *locations))
    
    def is_orphan(self) -> bool:
        """ Checks if the node is an orphan, that doesn't have any inbonding
        or outbonding edge. """
        return (not self.in_edges) and (not self.out_edges)
    
    def get_reference_next(self) -> Node:
        """ Get the next node of which the edge is reference (not variant
        or cleavage). """
        for edge in self.out_edges:
            if edge.type in ['reference', 'variant_end']:
                return edge.out_node
        return None
    
    def get_reference_prev(self) -> Node:
        """ Get the previous node of which the edge is reference (not variant
        or cleavage) """
        for edge in self.in_edges:
            if edge.type == 'reference':
                return edge.out_node
        return None

class Cursor():
    """ Defines the cursor when walking through a TranscriptVariantGraph
    object.

    Attributes:
        node (Node): The current node pointing at.
        seq (DNASeqRecordWithCoordinates): The current sequence.
        position_to_orf (str): Relative postion to the ORF. 'up' for upstream,
            meaning it's before the start codon, 'mid' for midstream meaning
            it's between the start and stop codon, and 'down' for downstream
            meaning the stop codon is in the sequence.
        search_orf (bool): Whether to search ORF. Only applied to sequances
            with a frameshifting in the upstream.
        cleave_sites (List[int]): The enzymatic cleave sites of the current
            sequence. Starts from 0.
        variants (List[VEPVariantRecord]): The variants that the current
            sequence contains.
        frameshifting (List[VEPVariantRecord]): The frameshifting mutations
            to the upstream.
    """
    def __init__(self, node:Node, seq:DNASeqRecordWithCoordinates,
            in_edge:Edge=None, position_to_orf:str='up', search_orf:bool=False,
            cleave_sites:List[int]=None, variants:list[VEPVariantRecord]=None,
            frameshifting:list[VEPVariantRecord]=None):
        """ Constructor for Cursor. """
        if cleave_sites is None:
            cleave_sites = []
        if variants is None:
            variants = []
        if frameshifting is None:
            frameshifting = []
        self.node = node
        self.in_edge = in_edge
        self.seq = seq
        self.position_to_orf = position_to_orf
        self.search_orf = search_orf
        self.cleave_sites = cleave_sites
        self.variants = variants
        self.frameshifting = frameshifting

    def __getitem__(self, index) -> Cursor:
        """ Get item. Only the sequence is subsetted. The variants and cleave
        sites are also removed if they are not in the range any more. """
        start,stop,step = index.indices(len(self.seq))
        seq = self.seq[start:stop]

        cleave_sites = []
        for cleave_site in self.cleave_sites:
            site = cleave_site - start
            if site > 0:
                cleave_sites.append(site)

        # only keep  the variants that are still in the range.
        # NOTE(CZ): Because neucleotides caused by mutations are not recorded
        # in the MatchedLocation. And if the sequence starts with a SNP, then
        # the variant location is actually before the first location of seq.
        if len(seq.locations) == 0:
            raise ValueError('sequence has no locations')
        variants = []
        for variant in self.variants:
            if seq.contains_variant(variant):
                variants.append(variant)
        
        return self.__class__(
            node=self.node,
            seq=seq,
            position_to_orf=self.position_to_orf,
            search_orf=self.search_orf,
            cleave_sites=cleave_sites,
            variants=variants,
            frameshifting=self.frameshifting
        )


    def join(self, node:Node, in_edge:Edge, start:int=0,
            position_to_orf:str=None, search_orf:bool=None,
            check_linked:bool=False) -> Cursor:
        """ Join the current cursor with the next. The sequence is also
        combined, as well as variants and cleave sites. """
        if check_linked:
            linked = False
            for edge in self.node.out_edges:
                if node is edge.out_node:
                    linked = True
            if not linked:
                raise ValueError('node is not out bonded.')
        
        tail = self[start:]
        if position_to_orf is None:
            position_to_orf = tail.position_to_orf
        if search_orf is None:
            search_orf = tail.search_orf
        
        variants = tail.variants
        frameshifting = tail.frameshifting
        if node.variant:
            variants.append(node.variant)
            if node.variant.is_frameshifting():
                frameshifting.append()

        return self.__class__(
            node=node,
            seq=tail.seq + node.seq,
            position_to_orf=position_to_orf,
            search_orf=search_orf,
            cleave_sites=tail.cleave_sites,
            variants=variants,
            frameshifting=frameshifting
        )
    
    def search_more_cleave_sites(self, rule:str, exception:str=None):
        """ Search for cleave sites after the last site. """
        start = 0 if len(self.cleave_sites) == 0 else self.cleave_sites[-1]
        end = len(self.seq)
        end = end - (end - start) % 3
        if end <= start:
            return
        self.cleave_sites += self.seq.find_all_enzymatic_cleave_sites(
            rule=rule,
            exception=exception,
            start=start,
            end=end
        )
    
    def has_any_variants(self) -> bool:
        """ Returns if the current cursor node has any variants """
        return len(self.variants) > 0 or len(self.frameshifting) > 0

class TranscriptVariantGraph():
    """
    Defines the DAG(directed acyclic graph) data structure for the algorithm to
    find variant peptides raised by genomic variants.

    The algorithm has two main steps:

    1. Create the DAG with transcript variant records. For each variant, it
    first breaks the transcript into three nodes, as 'prior', 'between', and
    'post' to the variant location, and linked with edges that labeled as
    "reference". The variant is then linked to the pior and post node with
    edges labeled as 'variant_start' and 'variatn_end'.

    2. Walk through the graph, find enzymatic sites, translate to peptides,
    and add to the final results if it contains a variant and does not have
    the same sequence as any canonical peptides from the reference transcript
    sequence.

    Attributes:
        nodes (Set[Node]): Set of nodes. Each node ha
        edges (Set[Edge]): Set of edges. Each 
        seq (DNASeqRecordWithCoordinates): The original sequence of the 
            transcript (reference).
        transcript_id (str)
        variants List[VEPVariantRecord]: All variant records.
    """
    def __init__(self, seq:DNASeqRecordWithCoordinates, transcript_id:str):
        """ Constructor to create a TranscriptVariantGraph object.
        
        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the 
                transcript (reference).
            transcript_id (str)
        """
        self.V = 1
        if seq.locations == []:
            seq.locations = [MatchedLocation(
                query=FeatureLocation(start=0, end=len(seq)),
                end=FeatureLocation(start=0, end=len(seq))
            )]
        self.seq = seq
        node = Node(seq)
        self.root = node
        self.nodes = set([node])
        self.edges = set()
        self.transcript_id = transcript_id
        self.variants = []
    
    def variants_are_adjacent(self, edge1:Edge, edge2:Edge) -> bool:
        """ Deprecated If two edges are both variants. """
        if edge1 is None:
            return False
        if edge2 is None:
            return False
        return edge1.type == 'variant_start' and edge2.type == 'variant_start'
    
    def remove_edge(self, edge:Edge) -> None:
        """ Removes an edge from the graph.
        
        The edge will be removed from
        the outbond edge list of its inbond node, and from the inbond edge
        list of its outbond node. The inbond or outbond node is removed if it
        became an orphan (do not have any edge connected) after the removal
        of the edge.
        
        Args:
            edge (Edge): The edge to be removed.
        """
        if edge in edge.in_node.out_edges:
            edge.in_node.out_edges.remove(edge)
        if edge in edge.out_node.in_edges:
            edge.out_node.in_edges.remove(edge)
        if edge.in_node.is_orphan():
            self.nodes.remove(edge.in_node)
        if edge.out_node.is_orphan():
            self.nodes.remove(edge.out_node)
        if edge in self.edges:
            self.edges.remove(edge)
    
    def add_edge(self, in_node:Node, out_node:Node, type:str) -> Edge:
        """ Add an edge between two nodes.
        
        An edge is created abd added to the the outbond edge list of the inbond
        node, as well as the intbond edge list of the output node, then added
        to the edge list of the object.

        Args:
            in_node (Node): the inbond node of the edge to add.
            out_node (Node): the outbond node of the edge to add.
            type (str): the type of the edge to add. Can be either 'reference',
                'variant_start', and 'variant_end'.
        
        Returns:
            The new added edge.
        """
        edge = Edge(in_node, out_node, type)
        self.edges.add(edge)
        in_node.out_edges.add(edge)
        out_node.in_edges.add(edge)
        self.nodes.add(in_node)
        self.nodes.add(out_node)
        return edge

    def splice(self, node:Node, i:int, type:str) -> Tuple[Node, Node]:
        """ Splice a node into two separate nodes.

        The sequence of a given node is spliced at the given position. Two
        nodes with smaller sequences are then created with a edge pointing
        from the upstream to the downstream node. Any node connected two the 
        original node is transfered two the new nodes.

        Args:
            node (Node): The node to be spliced.
            i (int): The position of sequence to splice.
            type (str): The type 
            type (str): the type of the edge to add bewen the two new nodes.
                Can be either 'reference', 'variant_start', and 'variant_end'.
        
        Returns:
            The two new nodes.
        """
        left_node = Node(node.seq[:i])
        right_node = Node(node.seq[i:])

        while node.in_edges:
            edge = node.in_edges.pop()
            self.add_edge(edge.in_node, left_node, type=edge.type)
            self.remove_edge(edge)
        
        while node.out_edges:
            edge = node.out_edges.pop()
            self.add_edge(right_node, edge.out_node, type=edge.type)
            self.remove_edge(edge)
        
        self.add_edge(left_node, right_node, type=type)

        if node in self.nodes:
            self.nodes.remove(node)
        
        self.V += 1
        return left_node, right_node
    
    
    def apply_variant(self, node:Node, variant:VEPVariantRecord) -> Node:
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
            seq=variant.variant.alt,
            locations=[],
            orf=None
        )
        var_node = Node(seq = seq, variant=variant)
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
        if variant_end <= node_end:
            index = body.seq.get_query_index(variant_end)
            mid, tail = self.splice(body, index, 'reference')
            self.add_edge(var_node, tail, 'variant_end')
        else:
            # This is the case that the range of this variant is larger than
            # the node, such as deletion.
            cur = body
            while cur.seq.locations[-1].ref.end < variant_end:
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

            if cur.seq.locations[-1].ref.end < variant.variant.location.start:
                try:
                    cur = cur.get_reference_next()
                except StopIteration:
                    break
                continue

            cur = self.apply_variant(cur, variant)
            if len(cur.in_edges) == 0:
                self.root = cur
            variant = next(variant_iter, None)
    
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
        ref_orf = self.root.seq.orf[0]
        protein = self.seq[ref_orf.start:ref_orf.end].translate()

        try:
            carnonical_peptides = protein.enzymatic_cleave(
                rule=rule,
                miscleavage=miscleavage,
                exception=exception,
                min_mw=min_mw
            )
        except ValueError:
            raise
        carnonical_peptides = set([peptide.seq for peptide in \
            carnonical_peptides])
        # DNASeqRecordWithCoordinates only hashes the sequence, so here we are
        # only keeping the unique sequences of the carnonical peptides,
        # regardless of the location.
        carnonical_peptides = set(carnonical_peptides)
        incarnonical_peptides = set()

        queue = deque()
        queue.appendleft(Cursor(self.root, self.root.seq))
        visited = set()

        while queue:
            cursor:Cursor = queue.pop()

            if cursor.position_to_orf == 'up':
                # when it passes the ref_orf start site, and found that there
                # is a start codon lost mutation.
                if cursor.search_orf:
                    orf_start = cursor.seq.find_orf_first()
                    for edge in cursor.node.out_edges:
                        node = edge.out_node
                        if orf_start == -1:
                            # It is truncated because the first part is already
                            # check that there isn't a start codon.
                            new_cursor = cursor.join(node, edge, -2, 'up', True)
                        else:
                            new_cursor = cursor.join(
                                node, edge, orf_start, 'mid', True)
                        queue.appendleft(new_cursor)
                        continue
                # when it is before the ref_orf start
                if cursor.seq.locations[-1].ref.end - 3 < ref_orf.start:
                    for edge in cursor.node.out_edges:
                        new_cursor = cursor.join(edge.out_node, edge, -2)
                        queue.appendleft(new_cursor)
                    continue
                # when it is after the ref_orf start, we then need to check if
                # start codon is indeed there
                for edge in cursor.node.out_edges:
                    index = cursor.seq.get_query_index(ref_orf.start)
                    new_cursor = cursor.join(edge.out_node, edge, index, 'mid')
                    if new_cursor.seq.seq[:3] != 'ATG':
                        new_cursor.search_orf = True
                        new_cursor.position_to_orf = 'up'
                    queue.appendleft(new_cursor)
                continue
            
            # If a frameshifting mutation occured after start codon, the
            # stop codon needs to be searched.
            if len(cursor.frameshifting) > 0 \
                    and not cursor.search_orf \
                    and cursor.position_to_orf != 'down':
                cursor.search_orf = True
            
            # In case of frameshifting mutation, that the ORF needs to be
            # searched
            if cursor.search_orf:
                stop_codon = seq.find_stop_codon()
                if stop_codon > -1 :
                    new_cursor = cursor[:stop_codon]
                    new_cursor.search_orf = False
                    new_cursor.position_to_orf = 'down'
                    queue.appendleft(new_cursor)
                    continue
            
            if cursor.position_to_orf == 'mid' and \
                    len(cursor.node.seq.locations) > 0 and \
                    cursor.node.seq.locations[-1].ref.end > ref_orf.end:
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
                for edge in cursor.node.out_edges:
                    new_cursor = cursor.join(edge.out_node, edge)
                    new_cursor.search_more_cleave_sites(rule=rule,
                            exception=exception)
                    queue.appendleft(new_cursor)
                continue
            
            if not cursor.has_any_variants():
                for edge in cursor.node.out_edges:
                    new_cursor = cursor.join(edge.out_node, edge,
                        cursor.cleave_sites[0])
                    queue.appendleft(new_cursor)
                continue
            # cleave_iter = cursor.seq.iter_enzymatic_cleave_sites(
            #     rule=rule, exception=exception)
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
                    break
                seq = cursor.seq[:cleave_site]
                
                # check if the sequence is already visited
                visited_len_before = len(visited)
                seq_data_to_visit = str(seq.seq)
                visited.add(seq_data_to_visit)
                visited_len_after = len(visited)
                if visited_len_after == visited_len_before:
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
                
                # remove peptides without variants
                variants = []
                for variant in cursor.variants:
                    if seq.contains_variant(variant):
                        variants.append(variant)
                if len(variants) == 0:
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
