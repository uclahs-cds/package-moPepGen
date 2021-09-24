""" Module for transcript (DNA) variant graph """
from __future__ import annotations
from typing import List, Tuple, Set, Deque, Union, TYPE_CHECKING
from collections import deque
import copy
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen import dna, seqvar, svgraph
from moPepGen.svgraph.TVGCursor import TVGCursor


if TYPE_CHECKING:
    from moPepGen import gtf

class TranscriptVariantGraph():
    """ Defines the DAG data structure for the transcript and its variants.

    Attributes:
        nodes (Set[svgraph.TVGNode]): Set of nodes. Each node ha
        edges (Set[svgraph.TVGEdge]): Set of edges. Each
        seq (DNASeqRecordWithCoordinates): The original sequence of the
            transcript (reference).
        id (str)
        variants List[seqvar.VariantRecordWithCoordinate]: All variant
            records.
    """
    def __init__(self, seq:Union[dna.DNASeqRecordWithCoordinates,None],
            _id:str, cds_start_nf:bool=False, has_known_orf:bool=True,
            reading_frames:List[svgraph.TVGNode]=None):
        """ Constructor to create a TranscriptVariantGraph object.

        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the
                transcript (reference).
            transcript_id (str)
        """
        self.seq = seq
        if self.seq and not self.seq.locations:
            self.add_default_sequence_locations()
        node = svgraph.TVGNode(seq)
        self.root = node
        self.id = _id
        self.cds_start_nf = cds_start_nf
        self.has_known_orf = has_known_orf
        if reading_frames and len(reading_frames) != 3:
            raise ValueError('The length of reading_frames must be exactly 3.')
        self.reading_frames = reading_frames or [None, None, None]

    def add_default_sequence_locations(self):
        """ Add default sequence locations """
        self.seq.locations = [MatchedLocation(
            query=FeatureLocation(start=0, end=len(self.seq)),
            ref=FeatureLocation(start=0, end=len(self.seq))
        )]

    @staticmethod
    def variants_are_adjacent(edge1:svgraph.TVGEdge, edge2:svgraph.TVGEdge
            ) -> bool:
        """ [Deprecated] If two edges are both variants. """
        if edge1 is None:
            return False
        if edge2 is None:
            return False
        return edge1.type == 'variant_start' and edge2.type == 'variant_start'

    @staticmethod
    def remove_edge(edge:svgraph.TVGEdge) -> None:
        """ Removes an edge from the graph.

        The edge will be removed from
        the outbond edge list of its inbond node, and from the inbond edge
        list of its outbond node. The inbond or outbond node is removed if it
        became an orphan (do not have any edge connected) after the removal
        of the edge.

        Args:
            edge (svgraph.TVGEdge): The edge to be removed.
        """
        if edge in edge.in_node.out_edges:
            edge.in_node.out_edges.remove(edge)
        if edge in edge.out_node.in_edges:
            edge.out_node.in_edges.remove(edge)

    @staticmethod
    def add_edge(in_node:svgraph.TVGNode, out_node:svgraph.TVGNode,
            _type:str) -> svgraph.TVGEdge:
        """ Add an edge between two nodes.

        An edge is created abd added to the the outbond edge list of the inbond
        node, as well as the intbond edge list of the output node, then added
        to the edge list of the object.

        Args:
            in_node (svgraph.TVGNode): the inbond node of the edge to add.
            out_node (svgraph.TVGNode): the outbond node of the edge to add.
            type (str): the type of the edge to add. Can be either 'reference',
                'variant_start', and 'variant_end'.

        Returns:
            The new added edge.
        """
        edge = svgraph.TVGEdge(in_node, out_node, _type)
        in_node.out_edges.add(edge)
        out_node.in_edges.add(edge)
        return edge

    def remove_node(self, node:svgraph.TVGNode):
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
        new_root = svgraph.TVGNode(None)
        self.root = new_root
        self.add_edge(self.root, original_root, 'reference')

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
            for edge in cur.out_edges:
                if edge.out_node in visited:
                    continue
                if edge.out_node is None:
                    continue
                queue.appendleft(edge.out_node)
        return k

    @staticmethod
    def get_orf_index(node:svgraph.TVGNode, i:int=0) -> int:
        """ Find the ORF index of a given node at given position of its
        sequence """
        return node.get_orf_start(i) % 3

    def update_frameshifts(self, node:svgraph.TVGNode) -> None:
        """ Update the frameshifts slot """
        if node is not self.root:
            if len(node.in_edges) == 1:
                upstream = next(iter(node.in_edges)).in_node
            else:
                upstream = node.get_reference_prev()
            if upstream:
                node.frameshifts = copy.copy(upstream.frameshifts)

        for variant in node.variants:
            if variant.variant.is_frameshifting():
                node.frameshifts.add(variant.variant)

    def splice(self, node:svgraph.TVGNode, i:int, _type:str
            ) -> Tuple[svgraph.TVGNode, svgraph.TVGNode]:
        """ Splice a node into two separate nodes.

        The sequence of a given node is spliced at the given position. Two
        nodes with smaller sequences are then created with a edge pointing
        from the upstream to the downstream node. Any node connected two the
        original node is transfered two the new nodes.

        Args:
            node (svgraph.TVGNode): The node to be spliced.
            i (int): The position of sequence to splice.
            type (str): The type
            type (str): the type of the edge to add bewen the two new nodes.
                Can be either 'reference', 'variant_start', and 'variant_end'.

        Returns:
            The two new nodes.
        """
        left_node = svgraph.TVGNode(node.seq[:i])
        right_node = svgraph.TVGNode(node.seq[i:])

        for variant in node.variants:
            if variant.location.start < i:
                left_node.variants.append(variant)
            if variant.location.end > i:
                right_node.variants.append(variant)

        for frameshift in node.frameshifts:
            left_node.frameshifts.add(frameshift)
            right_node.frameshifts.add(frameshift)

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


    def apply_variant(self, node:svgraph.TVGNode, variant:seqvar.VariantRecord
            ) -> svgraph.TVGNode:
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
            input node is returned.
        """
        variant_start = variant.location.start
        variant_end = variant.location.end
        node_start = node.seq.locations[0].ref.start
        node_end = node.seq.locations[-1].ref.end
        if variant_start < node_start or variant_start > node_end:
            raise ValueError('Variant out of range')
        if variant.type == 'Deletion':
            seq = self.seq[variant.location.start:variant.location.start+1]
        else:
            seq = dna.DNASeqRecordWithCoordinates(
                seq=Seq(variant.alt),
                locations=[],
                orf=None
            )
        variant_with_coordinates = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=0, end=len(seq))
        )
        var_node = svgraph.TVGNode(
            seq=seq,
            variants=[variant_with_coordinates]
        )
        head = None
        mid = None
        # variant start
        if variant_start == node_start:
            # No need to splice, just add a edge to the prev node This is the
            # case that another variant happend at the same location
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

    def apply_fusion(self, node:svgraph.TVGNode, variant:seqvar.VariantRecord,
            variant_pool:seqvar.VariantRecordPool, genome:dna.DNASeqDict,
            anno:gtf.GenomicAnnotation) -> svgraph.TVGNode:
        """ Apply a fusion variant, by creating a subgraph of the donor
        transcript and merge at the breakpoint position.

        Note that all peptides from the donor transcripts are counted as
        variant peptide at this stage, and are filtered out when comparing to
        the global canonical peptide pool. This can be improved.

        Args:
            node (svgraph.TVGNode): The node where the fusion variant should be
                added to.
            variant (seqvar.VariantRecord): The fusion variant.
            donor_seq (dna.DNASeqRecordWithCoordinates): The donor transcript's
                sequence.
            donor_variants (List[seqvar.VariantRecord]): Variants that are
                associated with the donor transcript. Variants before the
                breakpoint won't be applied.
        """
        accepter_gene_id = variant.attrs['ACCEPTER_GENE_ID']
        accepter_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
        accepter_tx_model = anno.transcripts[accepter_tx_id]
        accepter_chrom = accepter_tx_model.transcript.location.seqname
        accepter_tx_seq = accepter_tx_model.get_transcript_sequence(genome[accepter_chrom])
        breakpoint_gene = variant.get_accepter_position()
        breakpoint_tx = anno.coordinate_gene_to_transcript(breakpoint_gene,
            accepter_gene_id, accepter_tx_id)

        accepter_variant_records = variant_pool.filter_variants(
            accepter_gene_id, anno, genome, exclude_type=['Fusion'],
            start=breakpoint_tx, return_coord='transcript'
        )

        branch = TranscriptVariantGraph(accepter_tx_seq, self.id)
        branch.root.seq = branch.root.seq[breakpoint_tx:]
        branch.add_null_root()
        branch.create_variant_graph(accepter_variant_records)
        var_node = branch.root.get_reference_next()
        while var_node.in_edges:
            edge = var_node.in_edges.pop()
            branch.remove_edge(edge)
        var_node.branch = True
        var = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=0, end=0)
        )
        var_node.variants.append(var)
        var_node.frameshifts.add(variant)
        variant_start = variant.location.start
        node_start = node.seq.locations[0].ref.start

        if variant_start < node_start:
            raise ValueError('Variant out of range')

        if variant_start == node_start:
            prev = node.get_reference_prev()
            self.add_edge(prev, var_node, 'variant_start')
            return node

        index = node.seq.get_query_index(variant_start)
        head, _ = self.splice(node, index, 'reference')
        self.add_edge(head, var_node, 'variant_start')
        return head

    def _apply_insertion(self, node:svgraph.TVGNode,
            var:seqvar.VariantRecordWithCoordinate,
            seq:dna.DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecord]
        ) -> svgraph.TVGNode:
        """ Apply insertion """
        branch = TranscriptVariantGraph(seq, self.id)
        branch.root.variants.append(var)
        if len(seq.seq) % 3 != 0:
            branch.root.frameshifts.add(var.variant)

        branch.add_null_root()
        branch.create_variant_graph(variants)
        var_head = branch.root.get_reference_next()
        var_tail = var_head
        next_node = var_tail.get_reference_next()
        while next_node:
            var_tail = next_node
            next_node = var_tail.get_reference_next()
        while var_head.in_edges:
            edge = var_head.in_edges.pop()
            branch.remove_edge(edge)
        var_head.branch = True

        var_start = var.variant.location.start
        var_end = var.variant.location.end
        node_start = node.seq.locations[0].ref.start
        node_end = node.seq.locations[-1].ref.end

        if var_start < node_start:
            raise ValueError('Variant out of range')

        if var_start == node_start:
            prev = node.get_reference_prev()
            self.add_edge(prev, var_head, 'variant_start')
            body = node
        else:
            index = node.seq.get_query_index(var_start)
            head, body = self.splice(node, index, 'reference')
            self.add_edge(head, var_head, 'variant_start')

        if var_end < node_end:
            index = body.seq.get_query_index(var_end)
            mid, tail = self.splice(body, index, 'reference')
            self.add_edge(var_tail, tail, 'variant_end')
        else:
            cur = body
            while cur.seq.locations[-1].ref.end <= var_end and cur.out_edges:
                cur = cur.get_reference_next()

            if cur.seq.locations[-1].ref.end > var_end:
                index = cur.seq.get_query_index(var_end)
                if index == 0:
                    self.add_edge(var_tail, cur, 'variant_end')
                else:
                    _, right = self.splice(cur, index, 'reference')
                    self.add_edge(var_tail, right, 'variant_end')

        if head is not None:
            return head
        if mid is not None:
            return mid
        return body

    def apply_insertion(self, node:svgraph.TVGNode,
            variant:seqvar.VariantRecord,
            variant_pool:seqvar.VariantRecordPool,
            genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation
            ) -> svgraph.TVGNode:
        """ Apply an insertion into the the TVG graph. """
        gene_id = variant.attrs['DONOR_GENE_ID']
        gene_model = anno.genes[gene_id]
        chrom = gene_model.chrom
        donor_start = variant.get_donor_start()
        donor_end = variant.get_donor_end()
        gene_seq = gene_model.get_gene_sequence(genome[chrom])
        insert_seq = gene_seq[donor_start:donor_end]
        exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion']
        insert_variants = variant_pool.filter_variants(gene_id, anno, genome,
            exclude_type, start=donor_start, end=donor_end)
        # add the reference sequence to the insertion sequence.
        ref_seq_location = MatchedLocation(
            query=FeatureLocation(start=0, end=1),
            ref=FeatureLocation(
                seqname=self.id,
                start=variant.location.start,
                end=variant.location.end
            )
        )
        ref_seq = dna.DNASeqRecordWithCoordinates(
            Seq(variant.ref),
            locations=[ref_seq_location]
        )
        insert_seq = ref_seq + insert_seq
        var = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=1, end=len(insert_seq.seq))
        )
        return self._apply_insertion(node, var, insert_seq, insert_variants)

    def apply_substitution(self, node:svgraph.TVGNode,
            variant:seqvar.VariantRecord,
            variant_pool:seqvar.VariantRecordPool,
            genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation
            ) -> svgraph.TVGNode:
        """ Apply a substitution variant into the graph """
        gene_id = variant.attrs['DONOR_GENE_ID']
        gene_model = anno.genes[gene_id]
        chrom = gene_model.chrom
        donor_start = variant.get_donor_start()
        donor_end = variant.get_donor_end()
        gene_seq = gene_model.get_gene_sequence(genome[chrom])
        sub_seq = gene_seq[donor_start:donor_end]
        exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion']
        sub_variants = variant_pool.filter_variants(gene_id, anno, genome,
            exclude_type, start=donor_start, end=donor_end)
        var = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=0, end=len(sub_seq.seq))
        )
        return self._apply_insertion(node, var, sub_seq, sub_variants)


    def create_variant_graph(self, variants:List[seqvar.VariantRecord]
            ) -> None:
        """ Create a variant graph.

        With a list of genomic variants, incorprate each variant into the
        graph.

        Args:
            variant [seqvar.VariantRecord]: The variant record.
        """
        variant_iter = iter(variants)
        variant = next(variant_iter, None)
        if self.root.seq:
            cur = self.root
        else:
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

    def copy_node(self, node:svgraph.TVGNode) -> svgraph.TVGNode:
        """ Create a copy of a node and connect to the same up and downstream

        Args:
            node (svgraph.TVGNode): The node to replicate.

        Return:
            The replicate of the node.
        """
        node_copy = svgraph.TVGNode(
            seq=copy.copy(node.seq),
            variants=copy.copy(node.variants),
            frameshifts=copy.copy(node.frameshifts),
            branch=node.branch
        )
        for edge in node.in_edges:
            self.add_edge(edge.in_node, node_copy, edge.type)
        for edge in node.out_edges:
            self.add_edge(node_copy, edge.out_node, edge.type)
        return node_copy


    def create_branch(self, node:svgraph.TVGNode) -> svgraph.TVGNode:
        """ Create a branch by making a deep copy of the specified node and
        all its children. The branch attribute of the copied node is set
        to True. """
        new_node:svgraph.TVGNode = node.deepcopy()
        edges = copy.copy(node.in_edges)
        while edges:
            edge = edges.pop()
            self.add_edge(edge.in_node, new_node, edge.type)
            self.remove_edge(edge)
        self.remove_node(node)
        new_node.branch = True
        return new_node

    def merge_with_outbonds(self, node:svgraph.TVGNode) -> List[svgraph.TVGNode]:
        """ For a given node, merge it with all its outbound nodes. """

        in_edges = copy.copy(node.in_edges)
        out_edges:Set[svgraph.TVGEdge] = copy.copy(node.out_edges)
        out_nodes = []
        while out_edges:
            edge = out_edges.pop()
            out_node = edge.out_node
            out_node.seq = node.seq + out_node.seq
            out_node.branch = node.branch
            out_node.orf = node.orf

            variants = copy.copy(node.variants)
            for variant in out_node.variants:
                variants.append(variant.shift(len(node.seq)))
            out_node.variants = variants

            for in_edge in in_edges:
                self.add_edge(in_edge.in_node, out_node, in_edge.type)

            self.update_frameshifts(out_node)

            out_nodes.append(out_node)

        self.remove_node(node)

        return out_nodes

    @staticmethod
    def find_additional_end_nodes(start:svgraph.TVGNode, end:svgraph.TVGNode
            ) -> List[svgraph.TVGNode]:
        """ Find additional end node of active regions, e.g. Fusion. Because
        Fusions are by default branches, so not processed during the variant
        node alignment. """
        end_nodes = []
        cur = start
        if cur.variants:
            cur = cur.get_reference_next()
        while cur is not end:
            for edge in cur.out_edges:
                out_node = edge.out_node
                if edge.type in ['reference', 'variant_end']:
                    next_node = edge.out_node
                for variant in out_node.variants:
                    if variant.variant.type == 'Fusion':
                        end_node = out_node.find_farthest_node_with_overlap()
                        end_nodes.append(end_node)
            cur = next_node
        return end_nodes

    def branch_or_skip_frameshifts(self, upstream:svgraph.TVGNode,
            branch_out_size:int=100, max_frameshift_dist:int=60,
            max_frameshift_num:int=3
            ) -> List[svgraph.TVGNode]:
        """ Create a branch for any frameshifting mutation within the active
        region. If the created branch also has frameshifting mutations within
        its active region, they will be further branched out, too. """
        dist = max_frameshift_dist
        num = max_frameshift_num
        start_node = upstream
        main_end_node = upstream.find_farthest_node_with_overlap()
        node_to_branch = start_node.next_node_to_branch_out(
            to_node=main_end_node,
            branch_out_size=branch_out_size
        )
        end_nodes = []
        while node_to_branch:
            self.update_frameshifts(node_to_branch)
            if node_to_branch.should_skip_frameshift(dist, num):
                self.remove_node(node_to_branch)
            else:
                branched_node = self.create_branch(node_to_branch)
                self.update_frameshifts(branched_node)
                branch_ends = self.branch_or_skip_frameshifts(
                    branched_node, branch_out_size, dist, num)
                end_nodes.extend(branch_ends)

            main_end_node = upstream.find_farthest_node_with_overlap()
            node_to_branch = start_node.next_node_to_branch_out(
                to_node=main_end_node,
                branch_out_size=branch_out_size
            )
        end_nodes.append(main_end_node)
        additional_end_nodes = self.find_additional_end_nodes(start_node,
            main_end_node)
        end_nodes.extend(additional_end_nodes)
        return end_nodes

    def align_variants(self, node:svgraph.TVGNode, branch_out_size:int=100,
            branch_out_frameshifting:bool=True, max_frameshift_dist:int=60,
            max_frameshift_num:int=3) -> svgraph.TVGNode:
        r""" Aligns all variants at that overlaps to the same start and end
        position. Frameshifting mutations will be brached out

        For example:

                 T--                     T-G-CCCT
                /   \                   /
            ATGG-TCT-G-CCCT   ->    ATGG-TCTG-CCCT
                    \ /                 \    /
                     A                   TCTA

        Args:
            node (svgraph.TVGNode): The node of which the outbound nodes will
                be aligned.
            branch_out_size (int): The size limit that if a variant is larger
                that it, it will branch out even if it's not a frameshifting
                mutation.
        Returns:
            The original input node.
        """
        start_node = node
        if branch_out_frameshifting:
            end_nodes = self.branch_or_skip_frameshifts(node,
                branch_out_size, max_frameshift_dist, max_frameshift_num)
        else:
            end_nodes = [node.find_farthest_node_with_overlap()]

        new_nodes = set()
        queue = deque()
        trash = set()

        # start by removing the out edges from the start node.
        while start_node.out_edges:
            edge = start_node.out_edges.pop()
            out_node:svgraph.TVGNode = edge.out_node
            self.remove_edge(edge)
            queue.appendleft(out_node)
            new_nodes.add(out_node)

        while queue:
            cur:svgraph.TVGNode = queue.pop()
            if cur in end_nodes or not cur.out_edges:
                continue

            # because the new node will be reconstructed and added back
            new_nodes.remove(cur)

            out_edges = copy.copy(cur.out_edges)
            for out_edge in out_edges:
                out_node:svgraph.TVGNode = out_edge.out_node

                if out_node in end_nodes:
                    new_node = svgraph.TVGNode(
                        seq=cur.seq,
                        variants=copy.copy(cur.variants),
                        frameshifts=copy.copy(cur.frameshifts),
                        branch=cur.branch
                    )
                    for edge in cur.out_edges:
                        self.add_edge(new_node, edge.out_node, _type=edge.type)
                    new_nodes.add(new_node)
                    continue

                trash.add(out_node)

                # If the next node has any variants, shift the location by the
                # length of the cur node's sequence.
                new_variants = copy.copy(cur.variants)
                for variant in out_node.variants:
                    new_variants.append(variant.shift(len(cur.seq)))

                # create new node with the combined sequence
                frameshifts = copy.copy(cur.frameshifts)
                frameshifts.update(out_node.frameshifts)
                new_node = svgraph.TVGNode(
                    seq=cur.seq + out_node.seq,
                    variants=new_variants,
                    frameshifts=frameshifts,
                    branch=cur.branch
                )

                edges = copy.copy(out_node.out_edges)
                for edge in edges:
                    edge_type = 'variant_end' if new_node.variants \
                        else 'reference'
                    self.add_edge(new_node, edge.out_node, _type=edge_type)

                if out_node not in end_nodes:
                    queue.appendleft(new_node)

                new_nodes.add(new_node)
            # now remove the cur node from graph
            self.remove_node(cur)

        for trash_node in trash:
            self.remove_node(trash_node)

        # add now nodes to the graph
        for new_node in new_nodes:
            if new_node.variants:
                in_edge_type = 'variant_start'
            else:
                in_edge_type = 'reference'
            self.add_edge(start_node, new_node, in_edge_type)

        return start_node

    def expand_alignments(self, node:svgraph.TVGNode) -> List[svgraph.TVGNode]:
        r""" Expand the aligned variants into the range of codons. For
        frameshifting mutations, a copy of each downstream node will be
        created and branched out. Exclusive nodes are merged.

        Sequence after stop codons are trimmed off. However for frameshifting
        mutations where no stop codon is found, no trimming will be done. So
        the trailing sequence may not be multiple of 3.

        For example:

             T-GCCCT                GTGCCCT
            /                      /
        ATGG-TCTG-CCCT    ->    ATG-GTCTGC-CCT
            \    /                 \      /
             TCTA                   GTCTAC

        Return:
            A list of svgraph.TVGNode. For single sibling, the node after
            expension is returned, and for siblings more than 1, the common
            downstream node is returned.
        """
        # number of NT to be carried over to the downstreams.
        if node.seq:
            left_index = len(node.seq) - len(node.seq) % 3
            left_over = node.truncate_right(left_index)
        else:
            left_over_seq = dna.DNASeqRecordWithCoordinates(Seq(''), [])
            left_over = svgraph.TVGNode(left_over_seq)

        end_nodes = []

        # The reason of grouping outbond nodes into siblings
        for sibling_nodes in self.group_outbond_siblings(node):

            for sibling in sibling_nodes:
                sibling.append_left(left_over)
                sibling.orf = node.orf

            if len(sibling_nodes) == 1:
                end_node = self.single_expand_to_codons(sibling_nodes[0])
            else:
                end_node = self.siblings_expand_to_codons(sibling_nodes)

            if end_node:
                end_node.orf = node.orf
                end_nodes.append(end_node)

        return end_nodes

    def truncate_at_stop_codon_if_found(self, node:svgraph.TVGNode):
        """ Look for stop codon in the given node. If a stop codon is found,
        the node is then truncated at the stop codon position, and the
        downstream nodes will be removed. """
        stop_codon_index = node.seq.find_stop_codon()
        if stop_codon_index > -1:
            node.truncate_right(stop_codon_index)
            while node.out_edges:
                edge = node.out_edges.pop()
                self.remove_edge(edge)

    def siblings_expand_to_codons(self, nodes:List[svgraph.TVGNode]
            ) -> svgraph.TVGNode:
        """ Expand multiple sibling nodes to codons """
        sibling_len = len(nodes[0].seq.seq)

        ref_node = nodes[0].get_reference_next()
        for node in nodes:
            if not node.is_inbond_of(ref_node):
                raise ValueError("Sibling alignments don't have the same "
                    "outbound node.")

        right_index = (3 - sibling_len % 3) % 3

        if len(ref_node.seq.seq) < right_index + 3:
            # shouldn't happen. Raising an  exception for now.
            if len(ref_node.out_edges) > 0:
                raise ValueError('Somthing went wrong with variant alignment.')
            for node in nodes:
                node.append_right(ref_node)
                if self.has_known_orf:
                    self.truncate_at_stop_codon_if_found(node)
            self.remove_node(ref_node)
            return None

        found_stop = []
        right_over = ref_node.truncate_left(right_index)
        for node in nodes:
            node.append_right(right_over)
            if self.has_known_orf:
                self.truncate_at_stop_codon_if_found(node)
            found_stop.append(not node.out_edges)
        if all(found_stop):
            return None
        return ref_node

    def single_expand_to_codons(self, node:svgraph.TVGNode) -> svgraph.TVGNode:
        """ Expand a single node to codons """
        # If the given node is caused by a Fusion, we don't need to do anything
        # Because the acceptor (downstream) node is actually a reference node
        # of the acceptor transcript.
        if any(variant.variant.type == 'Fusion' for variant in node.variants):
            return node
        cur = node
        while len(cur.out_edges) == 1:
            cur = self.merge_with_outbonds(cur)[0]
        if len(cur.seq.seq) > 3 or not cur.out_edges:
            return cur
        # This condition shouldn't happen. I don't see how, but raising an
        # exception for now.
        raise ValueError('Something went wrong with variant alignment.')


    @staticmethod
    def group_outbond_siblings(upstream:svgraph.TVGNode
            ) -> List[List[svgraph.TVGNode]]:
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
        out_nodes = {edge.out_node for edge in upstream.out_edges}
        groups:List[List[svgraph.TVGNode]] = []
        while out_nodes:
            out_node = out_nodes.pop()
            downstream = out_node.get_reference_next()
            group:List[svgraph.TVGNode] = [out_node]
            if downstream:
                for edge in downstream.in_edges:
                    if edge.in_node in out_nodes:
                        group.append(edge.in_node)
                        out_nodes.remove(edge.in_node)
            groups.append(group)
        return groups

    def prune_variants_or_branch_out(self, node:svgraph.TVGNode
            ) -> Tuple[svgraph.TVGNode, List[svgraph.TVGNode]]:
        """ For a given reference node, remove the outbound variants if they
        are before the start codon. For start lost variants, create a branch.
        """
        outbound_ref = node.get_reference_next()
        outbound_ref.orf = node.orf
        ref_start = outbound_ref.seq.locations[0].ref.start
        orf = self.seq.orf
        branches = []

        # Maybe improved for start retained mutations
        if orf.start > ref_start:
            edges = copy.copy(node.out_edges)
            for edge in edges:
                if edge.out_node is outbound_ref:
                    continue
                self.remove_node(edge.out_node)

        edges = copy.copy(node.out_edges)
        for edge in edges:
            if edge.out_node is outbound_ref:
                continue
            self.update_frameshifts(edge.out_node)
            branch = self.create_branch(edge.out_node)
            branches.append(branch)
        if node is not self.root:
            nodes = self.merge_with_outbonds(node)
            for merged_node in nodes:
                if not merged_node.variants:
                    return merged_node, branches
        return outbound_ref, branches

    def skip_nodes_or_branch_out(self, node:svgraph.TVGNode, is_root:bool=False
            ) -> Tuple[svgraph.TVGNode, List[svgraph.TVGNode]]:
        """ For a given reference node, remove it if no start codon is found.
        Its outbond nodes are also removed if no start codon is found on them.
        Branch out frameshifting variants.
        """
        branches = []
        upstream_edges = copy.copy(node.in_edges)
        downstream = node.find_farthest_node_with_overlap()
        if is_root:
            nodes = [edge.out_node for edge in node.out_edges]
        else:
            nodes = self.merge_with_outbonds(node)

        # this is when all variants at this positions are frameshifting, or
        # downstream is None, or there is a variant at the last nucleotide.
        if downstream in nodes or downstream is None:
            for anode in nodes:
                if any(v.variant.is_frameshifting() for v in anode.variants):
                    if not anode.branch:
                        anode = self.create_branch(anode)
                if anode.branch:
                    branches.append(anode)
            return downstream, branches

        for anode in nodes:
            seq = anode.seq + downstream.seq[:2]
            index = seq.find_start_codon()
            if index == -1:
                self.remove_node(anode)
            else:
                self.update_frameshifts(anode)
                anode = self.create_branch(anode)
                anode.truncate_left(index)
                anodes = self.merge_with_outbonds(anode)
                if len(anodes) > 1:
                    raise ValueError('Multiple output nodes found.')
                branches.append(anodes[0])
        for edge in upstream_edges:
            self.add_edge(edge.in_node, downstream, _type=edge.type)
        return downstream, branches

    def find_orf_unknown(self) -> None:
        r""" Find all possible start codons when the ORF is unknown. Each ORF
        becomes the start of a subgraph downstream directly to the root.

        This is used for untranslatable transcripts such as lncRNA.

                                 ATGCT
                                /
        GGATGG-G-TGCT    ->    ^-ATGG-G-TGCT
              \ /                    \ /
               A                      A
        """
        cur = TVGCursor(node=self.root, search_orf=False)
        queue:Deque[TVGCursor] = deque([cur])

        while queue:
            cur = queue.pop()
            if cur.node.is_reference() and \
                    all(frame for frame in self.reading_frames):
                if self.root.is_inbond_of(cur.node):
                    self.remove_node(cur.node)
                continue
            if cur.node is self.root:
                if len(cur.node.out_edges) == 1:
                    main = cur.node.get_reference_next()
                    new_cursor = TVGCursor(node=main)
                    queue.append(new_cursor)
                    continue
                self.align_variants(cur.node)
                main, branches = self.skip_nodes_or_branch_out(cur.node, True)
                new_cursor = TVGCursor(node=main)
                queue.appendleft(new_cursor)
                for branch in branches:
                    new_cursor = TVGCursor(branch, True)
                    queue.appendleft(new_cursor)
                continue

            node = cur.node
            index = node.seq.find_start_codon()

            if index == -1:
                if not node.out_edges:
                    self.remove_node(node)
                    continue
                node.truncate_left(i=len(node.seq.seq)-2)
                if len(node.out_edges) == 1:
                    nodes = self.merge_with_outbonds(node)
                    new_cursor = TVGCursor(nodes[0])
                    queue.append(new_cursor)
                    continue
                node = self.align_variants(node)
                main, branches = self.skip_nodes_or_branch_out(node)
                if main:
                    new_cursor = TVGCursor(node=main)
                    queue.appendleft(new_cursor)
                for branch in branches:
                    new_cursor = TVGCursor(branch, True)
                    queue.appendleft(new_cursor)
                continue

            orf_start = node.get_orf_start(index)
            orf_id = orf_start % 3

            if cur.node.is_reference():
                if self.reading_frames[orf_id]:
                    node.truncate_left(index + 3)
                    new_cursor = TVGCursor(node)
                    queue.appendleft(new_cursor)
                    continue

            orf = [orf_start, None]
            node.truncate_left(index)
            node_copy = self.copy_node(node)
            self.update_frameshifts(node_copy)
            node_copy = self.create_branch(node_copy)

            if cur.node.is_reference():
                self.reading_frames[orf_id] = node_copy

            node_copy.orf = orf

            node.truncate_left(3)
            new_cursor = TVGCursor(node)
            queue.appendleft(new_cursor)

    def find_orf_known(self) -> None:
        """ Find the start codon according to the known ORF. Any variants
        occured prior to the start codon will be ignored. For start lost
        mutations, the start codon will be searched.
        """
        cur = TVGCursor(node=self.root, search_orf=False)
        queue:Deque[TVGCursor] = deque([cur])

        while queue:
            cur = queue.pop()
            if cur.node is self.root:
                if len(cur.node.out_edges) == 1:
                    main = cur.node.get_reference_next()
                    new_cursor = TVGCursor(node=main, search_orf=False)
                    queue.append(new_cursor)
                    continue
                main,branches = self.prune_variants_or_branch_out(cur.node)
                new_cursor = TVGCursor(node=main, search_orf=False)
                queue.appendleft(new_cursor)
                for branch in branches:
                    search_orf = not self.cds_start_nf
                    new_cursor = TVGCursor(branch, search_orf)
                    queue.appendleft(new_cursor)
                continue

            if cur.search_orf:
                node = cur.node
                index = node.seq.find_start_codon()

                # for start lost mutation or non-translatable sequences, if
                # start codon is not found in the current cursor, merge the
                # last nucleotides with each of the out nodes, and research.
                # This may be improved by skipping visited nodes.
                if index == -1:
                    if not node.out_edges:
                        self.remove_node(node)
                        continue
                    node.truncate_left(i=len(node.seq.seq)-2)
                    if len(node.out_edges) == 1:
                        nodes = self.merge_with_outbonds(node)
                        new_cursor = TVGCursor(nodes[0], True)
                        queue.appendleft(new_cursor)
                        continue
                    node = self.align_variants(node)
                    main, branches = self.skip_nodes_or_branch_out(node)
                    new_cursor = TVGCursor(node=main)
                    queue.appendleft(new_cursor)
                    for branch in branches:
                        new_cursor = TVGCursor(branch, True)
                        queue.appendleft(new_cursor)
                else:
                    orf_start = node.get_orf_start(index)
                    node.truncate_left(index)
                    node.orf = [orf_start, None]
                continue

            # back to the reference track
            if self.cds_start_nf:
                orf_start = 0
                start_codon_index = 0
            else:
                orf_start = int(self.seq.orf.start)
                start_codon_index = cur.node.seq.get_query_index(orf_start)

            if start_codon_index != -1:
                cur_node = cur.node
                cur_node.truncate_left(start_codon_index)
                cur_node.orf = [orf_start, None]
                while len(cur_node.seq) <= 2:
                    cur_node, branchs = self.prune_variants_or_branch_out(cur_node)
                    for branch in branchs:
                        new_cursor = TVGCursor(branch, True)
                        queue.appendleft(new_cursor)
                continue

            if not cur.node.out_edges:
                self.remove_node(cur.node)
                continue

            cur.node.truncate_left(i=len(cur.node.seq.seq)-2)
            main, branchs = self.prune_variants_or_branch_out(cur.node)
            new_cursor = TVGCursor(node=main, search_orf=False)
            queue.append(new_cursor)
            for branch in branchs:
                new_cursor = TVGCursor(node=branch, search_orf=True)
                queue.appendleft(new_cursor)

    def find_all_orfs(self):
        """ Fild all ORFs """
        if self.has_known_orf:
            self.find_orf_known()
        else:
            self.find_orf_unknown()

    def fit_into_codons(self, max_frameshift_dist:int=75,
            max_frameshift_num:int=3, branch_out_size:int=100) -> None:
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
        if self.root.seq is not None:
            self.add_null_root()
        queue:Deque[TVGCursor] = deque()
        for edge in self.root.out_edges:
            cur = TVGCursor(node=edge.out_node, search_orf=False)
            queue.append(cur)

        while queue:
            cur:TVGCursor = queue.pop()

            node = cur.node
            if not node.out_edges:
                continue

            # If the node has only one out node, we just need to merge then
            # together. This is the case of the reference node after a
            # frameshifting mutation is branched out.
            if len(node.out_edges) == 1:
                node = self.merge_with_outbonds(node)[0]
                if self.root.is_inbond_of(node):
                    orf_idx = node.orf[0] % 3
                    self.reading_frames[orf_idx] = node
                new_cursor = TVGCursor(node=node, search_orf=False)
                queue.appendleft(new_cursor)
                continue

            if len(node.out_edges) > 1:
                node = self.align_variants(node, branch_out_size,
                    max_frameshift_dist=max_frameshift_dist,
                    max_frameshift_num=max_frameshift_num)
            branches = self.expand_alignments(node)

            for branch in branches:
                new_cursor = TVGCursor(node=branch, search_orf=False)
                queue.appendleft(new_cursor)

    def translate(self) -> svgraph.PeptideVariantGraph:
        r""" Converts a DNA transcript variant graph into a peptide variant
        graph. A stop * is added to the end of all branches.

                GTG-CCCT              V-P-*
               /                     /
            ATG-GTCTGC-CCT    ->    M-VC-P-*
               \      /              \  /
                GTCTAC                VY
        """
        root = svgraph.PVGNode(None)
        pgraph = svgraph.PeptideVariantGraph(root, self.id, self.has_known_orf)

        queue = deque([(self.root, root)])
        visited = dict()

        while queue:
            dnode, pnode = queue.pop()
            if not dnode.out_edges:
                pgraph.add_stop(pnode)
                continue
            for edge in dnode.out_edges:
                out_node:svgraph.TVGNode = edge.out_node
                if out_node in visited:
                    out_pnode = visited[out_node]
                    pnode.add_out_edge(out_pnode)
                    continue

                orf = out_node.orf

                new_pnode = out_node.translate()
                new_pnode.orf = orf

                if self.has_known_orf:
                    stop_index = new_pnode.seq.seq.find('*')
                    if stop_index > -1:
                        new_pnode.truncate_right(stop_index)

                    pnode.add_out_edge(new_pnode)
                    visited[out_node] = new_pnode
                    if stop_index > -1:
                        pgraph.add_stop(new_pnode)
                    else:
                        if not out_node.out_edges:
                            new_pnode.truncated = True
                        queue.appendleft((out_node, new_pnode))
                else:
                    pnode.add_out_edge(new_pnode)
                    visited[out_node] = new_pnode
                    queue.appendleft((out_node, new_pnode))
        for i, dnode in enumerate(self.reading_frames):
            if dnode:
                pgraph.reading_frames[i] = visited[dnode]
        return pgraph
