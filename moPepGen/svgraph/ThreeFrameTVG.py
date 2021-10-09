""" Module for transcript (DNA) variant graph """
from __future__ import annotations
from typing import Dict, List, Tuple, Set, Deque, Union, TYPE_CHECKING
from collections import deque
import copy
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen import dna, seqvar, all_equal
from moPepGen.svgraph.TVGNode import TVGNode
from moPepGen.svgraph.TVGEdge import TVGEdge
from moPepGen.svgraph.PeptideVariantGraph import PeptideVariantGraph
from moPepGen.svgraph.PVGNode import PVGNode


if TYPE_CHECKING:
    from moPepGen import gtf

class ThreeFrameTVG():
    """ Defines the DAG data structure for the transcript and its variants.

    Attributes:
        seq (DNASeqRecordWithCoordinates): The original sequence of the
            transcript (reference).
        root (TVGNode)
        reading_frames (List[TVGNode])
        has_known_orf (bool)
        cds_start_nf (bool)
        id (str)
    """
    def __init__(self, seq:Union[dna.DNASeqRecordWithCoordinates,None],
            _id:str, root:TVGNode=None, reading_frames:List[TVGNode]=None,
            cds_start_nf:bool=False, has_known_orf:bool=None):
        """ Constructor to create a TranscriptVariantGraph object.

        Args:
            seq (DNASeqRecordWithCoordinates): The original sequence of the
                transcript (reference).
            transcript_id (str)
        """
        self.seq = seq
        if self.seq and not self.seq.locations:
            self.add_default_sequence_locations()
        self.root = root or TVGNode(None)
        self.id = _id
        self.reading_frames = reading_frames or [None, None, None]
        if reading_frames and len(reading_frames) != 3:
            raise ValueError('The length of reading_frames must be exactly 3.')
        self.cds_start_nf = cds_start_nf
        self.has_known_orf = has_known_orf or bool(seq.orf is not None)

    def add_default_sequence_locations(self):
        """ Add default sequence locations """
        self.seq.locations = [MatchedLocation(
            query=FeatureLocation(start=0, end=len(self.seq)),
            ref=FeatureLocation(start=0, end=len(self.seq))
        )]

    def init_three_frames(self, truncate_head:bool=True):
        """ Initiate the three reading frames as empty nodes.

        Args:
            truncated_head (bool): If true, the first x nucleotides are
                stripped off, to simulate how reading frames work. Defaults to
                True.
        """
        root0 = TVGNode(None, reading_frame_index=0)
        root1 = TVGNode(None, reading_frame_index=1)
        root2 = TVGNode(None, reading_frame_index=2)

        node0 = TVGNode(self.seq, reading_frame_index=0, subgraph_id=self.id)
        if truncate_head:
            node1 = TVGNode(self.seq[1:], reading_frame_index=1, subgraph_id=self.id)
            node2 = TVGNode(self.seq[2:], reading_frame_index=2, subgraph_id=self.id)
        else:
            node1 = TVGNode(self.seq, reading_frame_index=1, subgraph_id=self.id)
            node2 = TVGNode(self.seq, reading_frame_index=2, subgraph_id=self.id)

        self.add_edge(root0, node0, 'reference')
        self.add_edge(self.root, root0, 'reference')

        self.add_edge(root1, node1, 'reference')
        self.add_edge(self.root, root1, 'reference')

        self.add_edge(root2, node2, 'reference')
        self.add_edge(self.root, root2, 'reference')

        self.reading_frames = [root0, root1, root2]

    @staticmethod
    def remove_edge(edge:TVGEdge) -> None:
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
    def add_edge(in_node:TVGNode, out_node:TVGNode,
            _type:str) -> TVGEdge:
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
        edge = TVGEdge(in_node, out_node, _type)
        in_node.out_edges.add(edge)
        out_node.in_edges.add(edge)
        return edge

    def is_out_bond_to_any_root(self, node:TVGNode):
        """ Check if the node is out bond to any root of reading frame """
        for root in self.reading_frames:
            if root.is_inbond_of(node):
                return True
        return False

    def remove_node(self, node:TVGNode):
        """ Remove a node from its inboud and outbound node. """
        while node.in_edges:
            edge = node.in_edges.pop()
            self.remove_edge(edge)
        while node.out_edges:
            edge = node.out_edges.pop()
            self.remove_edge(edge)

    def get_known_reading_frame_index(self):
        """ get known reading frame index """
        return self.seq.orf.start % 3

    def should_add_singleton(self, node:TVGNode) -> bool:
        """ Checks if singleton frameshift should be added """
        return not self.has_known_orf or \
            node.reading_frame_index == self.get_known_reading_frame_index()

    def add_null_root(self):
        """ Adds a null node to the root. """
        original_root = self.root
        new_root = TVGNode(None)
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
    def get_orf_index(node:TVGNode, i:int=0) -> int:
        """ Find the ORF index of a given node at given position of its
        sequence """
        return node.get_orf_start(i) % 3

    def splice(self, node:TVGNode, i:int, _type:str
            ) -> Tuple[TVGNode, TVGNode]:
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
        left_node = TVGNode(
            node.seq[:i],
            subgraph_id=self.id,
            reading_frame_index=node.reading_frame_index,
            orf=node.orf
        )
        right_node = TVGNode(
            node.seq[i:],
            subgraph_id=self.id,
            reading_frame_index=node.reading_frame_index,
            orf=node.orf
        )

        for variant in node.variants:
            if variant.location.start < i:
                left_node.variants.append(variant)
            if variant.location.end > i:
                right_node.variants.append(variant)

        if node.frameshifts:
            left_node.frameshifts = copy.copy(node.frameshifts)
            right_node.frameshifts = copy.copy(node.frameshifts)
            for variant in node.variants:
                if variant.variant.is_frameshifting() and \
                        variant.location.start > len(left_node.seq.seq):
                    left_node.frameshifts.discard(variant)

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

        try:
            index = self.reading_frames.index(node)
            self.reading_frames[index] = left_node
        except ValueError:
            pass

        return left_node, right_node

    def apply_variant(self, source:TVGNode, target:TVGNode,
            variant:seqvar.VariantRecord) -> TVGNode:
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
        in_frame = (source is target)
        variant_start = variant.location.start
        variant_end = variant.location.end
        source_start = source.seq.locations[0].ref.start
        source_end = source.seq.locations[-1].ref.end
        target_start = target.seq.locations[0].ref.start
        target_end = target.seq.locations[-1].ref.end
        if variant_start < source_start or variant_start > source_end:
            raise ValueError('Variant out of source range of source')
        if variant_start < target_start or variant_start > target_end:
            raise ValueError('Variant out of source range of target')

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
        var_node = TVGNode(
            seq=seq,
            variants=[variant_with_coordinates],
            subgraph_id=self.id,
            reading_frame_index=source.reading_frame_index,
            orf=source.orf
        )
        returns = [None, None]
        # variant start
        if variant_start == source_start:
            # No need to splice, just add a edge to the prev node This is the
            # case that another variant happend at the same location
            prev = source.get_reference_prev()
            self.add_edge(prev, var_node, 'variant_start')
            returns[0] = source
        else:
            index = source.seq.get_query_index(variant_start)
            head, tail = self.splice(source, index, 'reference')
            self.add_edge(head, var_node, 'variant_start')
            returns[0] = head
            if in_frame:
                target = tail

        # varaint end
        if variant_end < target_end:
            index = target.seq.get_query_index(variant_end)
            head, tail = self.splice(target, index, 'reference')
            self.add_edge(var_node, tail, 'variant_end')
            returns[1] = head
        else:
            # This is the case that the range of this variant is larger than
            # the node, such as deletion.
            cur = target
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
            returns[1] = target

        if in_frame:
            returns[1] = returns[0]
        # returns the node with the first nucleotide of the input node.
        return returns

    def apply_fusion(self, cursors:List[TVGNode], variant:seqvar.VariantRecord,
            variant_pool:seqvar.VariantRecordPool, genome:dna.DNASeqDict,
            anno:gtf.GenomicAnnotation) -> List[TVGNode]:
        """ Apply a fusion variant, by creating a subgraph of the donor
        transcript and merge at the breakpoint position.

        Note that all peptides from the donor transcripts are counted as
        variant peptide at this stage, and are filtered out when comparing to
        the global canonical peptide pool. This can be improved.

        Args:
            node (TVGNode): The node where the fusion variant should be
                added to.
            variant (seqvar.VariantRecord): The fusion variant.
            donor_seq (dna.DNASeqRecordWithCoordinates): The donor transcript's
                sequence.
            donor_variants (List[seqvar.VariantRecord]): Variants that are
                associated with the donor transcript. Variants before the
                breakpoint won't be applied.
        """
        cursors = copy.copy(cursors)
        accepter_gene_id = variant.attrs['ACCEPTER_GENE_ID']
        accepter_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
        accepter_tx_model = anno.transcripts[accepter_tx_id]
        accepter_chrom = accepter_tx_model.transcript.location.seqname
        accepter_tx_seq = accepter_tx_model.get_transcript_sequence(
            genome[accepter_chrom])
        breakpoint_gene = variant.get_accepter_position()
        breakpoint_tx = anno.coordinate_gene_to_transcript(breakpoint_gene,
            accepter_gene_id, accepter_tx_id)

        accepter_variant_records = variant_pool.filter_variants(
            accepter_gene_id, anno, genome, exclude_type=['Fusion'],
            start=breakpoint_tx, return_coord='transcript'
        )

        branch = ThreeFrameTVG(accepter_tx_seq[breakpoint_tx:], self.id)
        branch.init_three_frames(truncate_head=False)
        for root in branch.reading_frames:
            list(root.out_edges)[0].out_node.subgraph_id = branch.id
        branch.create_variant_graph(accepter_variant_records, None, None, None)
        for edge in branch.root.out_edges:
            edge.out_node.variants.append(variant)

        for i in range(3):
            var_node = branch.reading_frames[i].get_reference_next()
            while var_node.in_edges:
                edge = var_node.in_edges.pop()
                branch.remove_edge(edge)
            var = seqvar.VariantRecordWithCoordinate(
                variant=variant,
                location=FeatureLocation(start=0, end=0)
            )
            var_node.variants.append(var)
            variant_start = variant.location.start

            node = cursors[i]
            node_start = node.seq.locations[0].ref.start

            if variant_start < node_start:
                raise ValueError('Variant out of range')

            if variant_start == node_start:
                prev = node.get_reference_prev()
                self.add_edge(prev, var_node, 'variant_start')
                continue

            index = node.seq.get_query_index(variant_start)
            node, _ = self.splice(node, index, 'reference')
            self.add_edge(node, var_node, 'variant_start')
            self.reading_frames[i] = node
            cursors[i] = node

        return cursors

    def _apply_insertion(self, cursors:List[TVGNode],
            var:seqvar.VariantRecordWithCoordinate,
            seq:dna.DNASeqRecordWithCoordinates,
            variants:List[seqvar.VariantRecord]
        ) -> List[TVGNode]:
        """ Apply insertion """
        cursors = copy.copy(cursors)
        branch = ThreeFrameTVG(seq, self.id)
        branch.init_three_frames(truncate_head=False)
        for edge in branch.root.out_edges:
            node = edge.out_node
            node.variants.append(var)
            is_frameshifting = False
            if var.variant.is_frameshifting():
                is_frameshifting = True
        branch.create_variant_graph(variants, None, None, None)

        for i in range(3):
            var_head = branch.reading_frames[i].get_reference_next()
            var_tail = var_head
            next_node = var_tail.get_reference_next()
            while next_node:
                var_tail = next_node
                next_node = var_tail.get_reference_next()
            while var_head.in_edges:
                edge = var_head.in_edges.pop()
                branch.remove_edge(edge)

            var_start = var.variant.location.start
            var_end = var.variant.location.end
            source = cursors[i]
            if is_frameshifting:
                frame_shifted = var.variant.frames_shifted()
                j = (frame_shifted + i) % 3
                target = cursors[j]
            else:
                target = source

            source_start = source.seq.locations[0].ref.start
            target_end = target.seq.locations[-1].ref.end

            if var_start < source_start:
                raise ValueError('Variant out of range')

            if var_start == source_start:
                prev = source.get_reference_prev()
                self.add_edge(prev, var_head, 'variant_start')
            else:
                index = source.seq.get_query_index(var_start)
                source, tail = self.splice(source, index, 'reference')
                self.add_edge(source, var_head, 'variant_start')
                if not is_frameshifting:
                    target = tail

            if var_end < target_end:
                index = target.seq.get_query_index(var_end)
                target, tail = self.splice(target, index, 'reference')
                self.add_edge(var_tail, tail, 'variant_end')
            else:
                cur = target
                while cur.seq.locations[-1].ref.end <= var_end and cur.out_edges:
                    cur = cur.get_reference_next()

                if cur.seq.locations[-1].ref.end > var_end:
                    index = cur.seq.get_query_index(var_end)
                    if index == 0:
                        self.add_edge(var_tail, cur, 'variant_end')
                    else:
                        _, right = self.splice(cur, index, 'reference')
                        self.add_edge(var_tail, right, 'variant_end')

            cursors[i] = source
            if is_frameshifting:
                cursors[j] = target

        return cursors

    def apply_insertion(self, cursors:List[TVGNode],
            variant:seqvar.VariantRecord,
            variant_pool:seqvar.VariantRecordPool,
            genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation
            ) -> List[TVGNode]:
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
        return self._apply_insertion(cursors, var, insert_seq, insert_variants)

    def apply_substitution(self, cursors:List[TVGNode],
            variant:seqvar.VariantRecord,
            variant_pool:seqvar.VariantRecordPool,
            genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation
            ) -> List[TVGNode]:
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
        return self._apply_insertion(cursors, var, sub_seq, sub_variants)


    def create_variant_graph(self, variants:List[seqvar.VariantRecord],
            variant_pool, genome, anno) -> None:
        """ Create a variant graph.

        With a list of genomic variants, incorprate each variant into the
        graph.

        Args:
            variant [seqvar.VariantRecord]: The variant record.
        """
        variant_iter = iter(variants)
        variant = next(variant_iter, None)
        cursors = copy.copy([x.get_reference_next() for x in self.reading_frames])

        if self.has_known_orf:
            start_codon = FeatureLocation(
                start=self.seq.orf.start,
                end=self.seq.orf.start + 3
            )

        while variant:
            if not any(cursors):
                break
            # skipping start lost mutations
            start_altering = variant.location.start < 3 or \
                (self.has_known_orf and variant.location.overlaps(start_codon))
            if start_altering:
                variant = next(variant_iter, None)
                continue

            if not any(any(y.is_inbond_of(x) for y in self.reading_frames) \
                    for x in cursors) and not \
                    all_equal([x.seq.locations[0].ref.start for x in cursors]):
                raise ValueError('Not all cursors are equal. Somthing messed up.')

            if any(c.seq.locations[0].ref.start > variant.location.start for c in cursors):
                variant = next(variant_iter, None)
                continue

            any_cursor_expired = False
            for i,cursor in enumerate(cursors):
                if cursor.seq.locations[-1].ref.end <= variant.location.start:
                    cursors[i] = cursors[i].get_reference_next()
                    any_cursor_expired = True
            if any_cursor_expired:
                continue

            if variant.type == 'Fusion':
                cursors = self.apply_fusion(cursors, variant, variant_pool, genome, anno)

            elif variant.type == 'Insertion':
                cursors = self.apply_insertion(cursors, variant, variant_pool, genome, anno)

            elif variant.type == 'Substitution':
                cursors = self.apply_substitution(cursors, variant, variant_pool, genome, anno)

            elif variant.is_frameshifting() and not start_altering:
                frames_shifted = variant.frames_shifted()
                for i in range(3):
                    j = (i + frames_shifted) % 3
                    cursors[i], cursors[j] = self.apply_variant(cursors[i],
                        cursors[j], variant)
            else:
                for i,_ in enumerate(cursors):
                    cursors[i], _ = self.apply_variant(cursors[i], cursors[i],
                        variant)
            variant = next(variant_iter, None)

    def copy_node(self, node:TVGNode) -> TVGNode:
        """ Create a copy of a node and connect to the same up and downstream

        Args:
            node (TVGNode): The node to replicate.

        Return:
            The replicate of the node.
        """
        node_copy = TVGNode(
            seq=copy.copy(node.seq),
            variants=copy.copy(node.variants),
            frameshifts=copy.copy(node.frameshifts),
            branch=node.branch,
            subgraph_id=node.subgraph_id
        )
        for edge in node.in_edges:
            self.add_edge(edge.in_node, node_copy, edge.type)
        for edge in node.out_edges:
            self.add_edge(node_copy, edge.out_node, edge.type)
        return node_copy

    def create_branch(self, node:TVGNode) -> TVGNode:
        """ Create a branch by making a deep copy of the specified node and
        all its children. The branch attribute of the copied node is set
        to True. """
        new_node:TVGNode = node.deepcopy()
        edges = copy.copy(node.in_edges)
        while edges:
            edge = edges.pop()
            self.add_edge(edge.in_node, new_node, edge.type)
            self.remove_edge(edge)
        self.remove_node(node)
        new_node.branch = True
        return new_node

    def merge_with_outbonds(self, node:TVGNode) -> List[TVGNode]:
        """ For a given node, merge it with all its outbound nodes. """

        in_edges = copy.copy(node.in_edges)
        out_edges:Set[TVGEdge] = copy.copy(node.out_edges)
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

            out_nodes.append(out_node)

        self.remove_node(node)

        return out_nodes

    @staticmethod
    def find_bridge_nodes_between(start:TVGNode, end=TVGNode
            ) -> Tuple[List[TVGNode], List[TVGNode]]:
        """ Find all bridge in and bridge out nodes between a start and a
        end node. """
        this_id = start.reading_frame_index
        bridge_in = []
        bridge_out = []
        visited = set()
        queue:Deque[TVGNode] = deque([e.out_node for e in start.out_edges])
        while queue:
            cur = queue.pop()
            len_before = len(visited)
            visited.add(cur)
            len_after = len(visited)
            if len_before == len_after:
                continue
            if any(e.out_node.reading_frame_index != this_id for e in cur.out_edges):
                bridge_out.append(cur)
                continue
            for e in cur.in_edges:
                if e.in_node.reading_frame_index != this_id:
                    bridge_in.append(e.in_node)
            if cur is not end:
                for e in cur.out_edges:
                    queue.appendleft(e.out_node)
        return bridge_in, bridge_out

    def align_variants(self, node:TVGNode) -> Tuple[TVGNode, TVGNode]:
        r""" Aligns all variants at that overlaps to the same start and end
        position. Frameshifting mutations will be brached out

        For example:

                 T--                     T-G-CCCT
                /   \                   /
            ATGG-TCT-G-CCCT   ->    ATGG-TCTG-CCCT
                    \ /                 \    /
                     A                   TCTA

        Args:
            node (TVGNode): The node of which the outbound nodes will
                be aligned.
            branch_out_size (int): The size limit that if a variant is larger
                that it, it will branch out even if it's not a frameshifting
                mutation.
        Returns:
            The original input node.
        """
        start_node = node
        end_node = node.find_farthest_node_with_overlap()
        end_nodes = [end_node]
        bridge_ins, bridge_outs = self.find_bridge_nodes_between(start_node, end_node)

        for bridge in bridge_outs:
            for e in bridge.out_edges:
                end_nodes.append(e.out_node)

        new_nodes = set()
        queue = deque()
        trash = set()
        bridge_map:Dict[TVGNode, TVGNode] = {x:x for x in bridge_ins}
        new_bridges:Set[TVGNode] = set()

        # start by removing the out edges from the start node.
        for edge in copy.copy(start_node.out_edges) :
            out_node = edge.out_node
            self.remove_edge(edge)
            queue.appendleft(out_node)

        for bridge in bridge_ins:
            queue.appendleft(bridge)

        while queue:
            cur:TVGNode = queue.pop()
            if cur in end_nodes or not cur.out_edges:
                continue

            for out_edge in copy.copy(cur.out_edges):
                out_node:TVGNode = out_edge.out_node

                if out_node in end_nodes:
                    new_node = TVGNode(
                        seq=cur.seq,
                        variants=copy.copy(cur.variants),
                        frameshifts=copy.copy(cur.frameshifts),
                        branch=cur.branch,
                        orf=cur.orf,
                        reading_frame_index=cur.reading_frame_index,
                        subgraph_id=cur.subgraph_id
                    )
                    trash.add(cur)
                    for edge in cur.out_edges:
                        self.add_edge(new_node, edge.out_node, _type=edge.type)
                    if cur not in bridge_map:
                        new_nodes.add(new_node)
                    else:
                        bridge_map[new_node] = bridge_map[cur]
                        new_bridges.add(new_node)
                    continue

                trash.add(out_node)

                # If the next node has any variants, shift the location by the
                # length of the cur node's sequence.
                new_variants = copy.copy(cur.variants)
                for variant in out_node.variants:
                    new_variants.append(variant.shift(len(cur.seq)))

                # create new node with the combined sequence
                frameshifts = copy.copy(cur.frameshifts)
                new_node = TVGNode(
                    seq=cur.seq + out_node.seq,
                    variants=new_variants,
                    frameshifts=frameshifts,
                    branch=cur.branch,
                    reading_frame_index=cur.reading_frame_index,
                    orf=cur.orf,
                    subgraph_id=cur.subgraph_id
                )

                for edge in copy.copy(out_node.out_edges):
                    edge_type = 'variant_end' if new_node.variants \
                        else 'reference'
                    self.add_edge(new_node, edge.out_node, _type=edge_type)

                if out_node not in end_nodes:
                    queue.appendleft(new_node)

                if cur in bridge_map:
                    bridge_map[new_node] = bridge_map[cur]
            # now remove the cur node from graph
            if cur not in bridge_map:
                self.remove_node(cur)

        # add now nodes to the graph
        for new_node in new_nodes:
            if new_node.variants:
                in_edge_type = 'variant_start'
            else:
                in_edge_type = 'reference'
            self.add_edge(start_node, new_node, in_edge_type)

        for bridge in new_bridges:
            original_bridge = bridge_map[bridge]
            for edge in original_bridge.in_edges:
                self.add_edge(edge.in_node, bridge, edge.type)
            trash.add(original_bridge)

        for trash_node in trash:
            self.remove_node(trash_node)

        return start_node, end_node

    @staticmethod
    def expand_alignments(start:TVGNode) -> TVGNode:
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
            A list of TVGNode. For single sibling, the node after
            expension is returned, and for siblings more than 1, the common
            downstream node is returned.
        """
        # number of NT to be carried over to the downstreams.
        if start.seq:
            left_index = len(start.seq) - len(start.seq) % 3
            left_over = start.truncate_right(left_index)
        else:
            left_over_seq = dna.DNASeqRecordWithCoordinates(Seq(''), [])
            left_over = TVGNode(
                left_over_seq,
                subgraph_id=start.subgraph_id
            )

        for out_edge in start.out_edges:
            out_node = out_edge.out_node
            out_node.append_left(left_over)


        ref_node = start.get_reference_next()
        if not ref_node.out_edges:
            return None

        end = ref_node.get_reference_next()
        right_index = (3 - len(ref_node.seq) % 3) % 3
        right_over = end.truncate_left(right_index)

        for in_edge in end.in_edges:
            in_node = in_edge.in_node
            in_node.append_right(right_over)

        return end

    def fit_into_codons(self) -> None:
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
        queue:Deque[TVGNode] = deque(self.reading_frames)
        while queue:
            cur = queue.pop()
            if not cur.out_edges:
                continue
            nxt = cur.get_reference_next()
            if len(cur.out_edges) == 1 and len(nxt.in_edges) == 1:
                if cur in self.reading_frames:
                    cur = cur.get_reference_next()
                else:
                    cur = self.merge_with_outbonds(cur)[0]
                queue.appendleft(cur)
                continue

            self.align_variants(cur)
            if cur.out_edges:
                node = self.expand_alignments(cur)
                queue.appendleft(node)

    def translate(self) -> PeptideVariantGraph:
        r""" Converts a DNA transcript variant graph into a peptide variant
        graph. A stop * is added to the end of all branches.

                GTG-CCCT              V-P-*
               /                     /
            ATG-GTCTGC-CCT    ->    M-VC-P-*
               \      /              \  /
                GTCTAC                VY
        """
        root = PVGNode(None, None)
        if self.has_known_orf:
            known_orf = [int(self.seq.orf.start), int(self.seq.orf.end)]
        else:
            known_orf = [None, None]

        pgraph = PeptideVariantGraph(
            root=root,
            _id=self.id,
            known_orf=known_orf
        )

        queue = deque([(dnode, root) for dnode in self.reading_frames])
        visited = dict()

        while queue:
            dnode, pnode = queue.pop()
            if not dnode.out_edges:
                pgraph.add_stop(pnode)
                continue
            for edge in dnode.out_edges:
                out_node:TVGNode = edge.out_node
                if out_node in visited:
                    out_pnode = visited[out_node]
                    pnode.add_out_edge(out_pnode)
                    continue

                if self.has_known_orf:
                    orf = known_orf
                else:
                    orf = out_node.orf
                new_pnode = out_node.translate()
                new_pnode.orf = orf
                pnode.add_out_edge(new_pnode)
                visited[out_node] = new_pnode
                queue.appendleft((out_node, new_pnode))
        for i, dnode in enumerate(self.reading_frames):
            dnode = dnode.get_reference_next()
            if dnode:
                pgraph.reading_frames[i] = visited[dnode]
        return pgraph