""" Module for transcript (DNA) variant graph """
from __future__ import annotations
from typing import Dict, List, Tuple, Set, Deque, Union, TYPE_CHECKING, Iterable,\
    Callable
from collections import deque
import copy
from Bio.Seq import Seq
from moPepGen import dna, seqvar, err, params, get_logger
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecordWithCoordinates
from moPepGen.dna.DNASeqRecord import DNASeqRecordWithCoordinates
from moPepGen.seqvar import VariantRecordWithCoordinate
from moPepGen.seqvar.VariantRecordPoolOnDisk import VariantRecordPoolOnDisk
from moPepGen.svgraph.TVGNode import TVGNode
from moPepGen.svgraph.TVGEdge import TVGEdge
from moPepGen.svgraph.PeptideVariantGraph import PeptideVariantGraph
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.svgraph.SubgraphTree import SubgraphTree
from moPepGen.svgraph.TVGNodeCollapser import TVGNodeCollapser


if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation
    from moPepGen.params import CleavageParams
    from moPepGen.dna import DNASeqDict
    from moPepGen.seqvar import VariantRecord

class ThreeFrameTVG():
    """ Defines the DAG data structure for the transcript and its variants.

    Attributes:
        seq (DNASeqRecordWithCoordinates): The original sequence of the
            transcript (reference).
        root (TVGNode): The root of the graph. The sequence of the root node
            must be None.
        reading_frames (List[TVGNode]): List of three null nodes. Each node is
            the root to the subgraph of the corresponding reading frame.
        has_known_orf (bool): whether
        cds_start_nf (bool)
        mrna_end_nf (bool)
        global_variant (VariantRecord)
        id (str)
    """
    def __init__(self, seq:Union[DNASeqRecordWithCoordinates,None],
            _id:str, root:TVGNode=None, reading_frames:List[TVGNode]=None,
            cds_start_nf:bool=False, has_known_orf:bool=None,
            mrna_end_nf:bool=False, global_variant:VariantRecord=None,
            coordinate_feature_type:str=None, coordinate_feature_id:str=None,
            subgraphs:SubgraphTree=None, hypermutated_region_warned:bool=False,
            cleavage_params:CleavageParams=None, gene_id:str=None,
            sect_variants:List[VariantRecordWithCoordinate]=None,
            max_adjacent_as_mnv:int=2):
        """ Constructor to create a TranscriptVariantGraph object. """
        self.seq = seq
        self.id = _id
        if self.seq and not self.seq.locations:
            self.add_default_sequence_locations()
        self.root = root or self.create_node(seq=None)
        self.reading_frames = reading_frames or [None, None, None]
        if reading_frames and len(reading_frames) != 3:
            raise ValueError('The length of reading_frames must be exactly 3.')
        self.cds_start_nf = cds_start_nf
        if has_known_orf is None:
            self.has_known_orf = bool(seq.orf is not None)
        else:
            self.has_known_orf = has_known_orf
        self.mrna_end_nf = mrna_end_nf
        self.global_variant = global_variant
        self.subgraphs = subgraphs or SubgraphTree()
        self.subgraphs.add_root(
            _id, variant=self.global_variant,
            feature_type=coordinate_feature_type,
            feature_id=coordinate_feature_id
        )
        self.hypermutated_region_warned = hypermutated_region_warned
        self.cleavage_params = cleavage_params or params.CleavageParams('trypsin')
        self.gene_id = gene_id
        self.sect_variants = sect_variants or []
        self.max_adjacent_as_mnv = max_adjacent_as_mnv

    def is_circ_rna(self) -> bool:
        """ If the graph is a circRNA """
        return self.global_variant and self.global_variant.is_circ_rna()

    def should_clip_trailing_nodes(self) -> bool:
        """ Checks whether the transcript (or circRNA) sequence trailing peptides
        should be avoided to be called as variant peptides. """
        return self.is_circ_rna() or self.mrna_end_nf

    def add_default_sequence_locations(self):
        """ Add default sequence locations """
        self.seq.locations = [MatchedLocation(
            query=FeatureLocation(start=0, end=len(self.seq)),
            ref=FeatureLocation(start=0, end=len(self.seq), seqname=self.id)
        )]

    def update_node_level(self, level:int):
        """ Update the level of all nodes """
        queue:Deque[TVGNode] = deque([self.root])
        visited:Set[TVGNode] = set()
        while queue:
            cur = queue.pop()
            if cur in visited:
                continue
            visited.add(cur)
            cur.level = level
            for edge in cur.out_edges:
                queue.appendleft(edge.out_node)

    def init_three_frames(self, truncate_head:bool=True):
        """ Initiate the three reading-frame graph.

        Args:
            truncated_head (bool): If true, the first x nucleotides are
                stripped off, to simulate how reading frames work. Defaults to
                True.
        """
        level = self.root.level
        root0 = TVGNode(None, reading_frame_index=0, level=level)
        root1 = TVGNode(None, reading_frame_index=1, level=level)
        root2 = TVGNode(None, reading_frame_index=2, level=level)

        node0 = TVGNode(
            seq=copy.deepcopy(self.seq), reading_frame_index=0, subgraph_id=self.id,
            global_variant=self.global_variant, level=level
        )
        if truncate_head:
            node1 = TVGNode(
                seq=self.seq[1:], reading_frame_index=1, subgraph_id=self.id,
                global_variant=self.global_variant, level=level
            )
            node2 = TVGNode(
                self.seq[2:], reading_frame_index=2, subgraph_id=self.id,
                global_variant=self.global_variant, level=level
            )
        else:
            node1 = TVGNode(
                copy.deepcopy(self.seq), reading_frame_index=1, subgraph_id=self.id,
                global_variant=self.global_variant, level=level
            )
            node2 = TVGNode(
                copy.deepcopy(self.seq), reading_frame_index=2, subgraph_id=self.id,
                global_variant=self.global_variant, level=level
            )

        for i,node in enumerate([node0, node1, node2]):
            for loc in node.seq.locations:
                loc.query.reading_frame_index = i

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

    @staticmethod
    def is_reference_edge(in_node:TVGNode, out_node:TVGNode):
        """ checks if this is a reference edge """
        if out_node.is_reference() and out_node.subgraph_id == in_node.get_last_subgraph_id():
            return True
        in_vars = {x.variant for x in in_node.variants}
        out_vars = {x.variant for x in out_node.variants}

        if in_vars == out_vars:
            return True

        if not in_vars and len(out_vars) == 1 and list(out_vars)[0].is_fusion():
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

    def get_known_reading_frame_index(self) -> int:
        """ get known reading frame index """
        return self.seq.orf.start % 3

    def should_add_singleton(self, node:TVGNode) -> bool:
        """ Checks if singleton frameshift should be added """
        return not self.has_known_orf or \
            node.reading_frame_index == self.get_known_reading_frame_index()

    def add_null_root(self):
        """ Adds a null node to the root. """
        original_root = self.root
        new_root = self.create_node(seq=None, level=original_root.level)
        self.root = new_root
        self.add_edge(self.root, original_root, 'reference')

    def add_stop_node(self, node:TVGNode):
        """ Add stop node after the given node """
        for edge in node.out_edges:
            out_node = edge.out_node
            if out_node.seq.seq == '' and \
                    out_node.reading_frame_index == node.reading_frame_index:
                return
        stop = TVGNode(
            DNASeqRecordWithCoordinates('', []),
            reading_frame_index=node.reading_frame_index
        )
        self.add_edge(node, stop, 'reference')

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

    def find_node(self, func:Callable) -> List[TVGNode]:
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
                for out_node in cur.get_out_nodes():
                    if not out_node in visited:
                        queue.append(out_node)
        return targets

    def gather_sect_variants(self, anno:GenomicAnnotation):
        """ Create selenocysteine trunction map """
        sect_variants:List[VariantRecordWithCoordinate] = []
        for sec in self.seq.selenocysteine:
            sect_var = seqvar.create_variant_sect(anno, self.id, sec.start)
            sect_loc = FeatureLocation(sec.start, sec.end)
            sect_variants.append(seqvar.VariantRecordWithCoordinate(
                location=sect_loc, variant=sect_var
            ))
        self.sect_variants = sect_variants

    @staticmethod
    def get_orf_index(node:TVGNode, i:int=0) -> int:
        """ Find the ORF index of a given node at given position of its
        sequence """
        return node.get_orf_start(i) % 3

    def create_node(self, seq:DNASeqRecordWithCoordinates,
            variants:List[VariantRecordWithCoordinate]=None,
            branch:bool=False, orf:List[int]=None, reading_frame_index:int=None,
            subgraph_id:str=None, level:int=0, global_variant:VariantRecord=None
            ) -> TVGNode:
        """ create a node """
        return TVGNode(
            seq=seq,
            variants=variants,
            branch=branch,
            orf=orf,
            reading_frame_index=reading_frame_index,
            subgraph_id=subgraph_id or self.id,
            level=level,
            global_variant=global_variant
        )

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
        left_node = node
        right_node = left_node.truncate_right(i)

        for edge in copy.copy(left_node.out_edges):
            self.add_edge(right_node, edge.out_node, edge.type)
            self.remove_edge(edge)

        self.add_edge(left_node, right_node, _type)

        if self.root is node:
            self.root = left_node

        try:
            index = self.reading_frames.index(node)
            self.reading_frames[index] = left_node
        except ValueError:
            pass

        return left_node, right_node

    def apply_variant(self, source:TVGNode, target:TVGNode,
            variant:VariantRecord) -> TVGNode:
        """ Apply a given variant to the graph.

        For a given genomic variant with the coordinates of the transcript,
        the reference transcript sequence is first spliced into three fragments
        (nodes) as 'prior', 'between', and 'post' to the variant location.
        Unless the pior or post are already spliced. A new node with the
        sequence of the alt sequence is created and linked to the 'prior' and
        'post' nodes.

        Args:
            node [node]: The node where the variant to be add.
            variant [VariantRecord]: The variant record.

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

        is_deletion = variant.type == 'Deletion'
        if is_deletion:
            seq = self.seq[variant.location.start:variant.location.start+1]
            for loc in seq.locations:
                loc.query.reading_frame_index = source.reading_frame_index
        else:
            seq = dna.DNASeqRecordWithCoordinates(
                seq=Seq(variant.alt),
                locations=[],
                orf=None
            )
        variant_with_coordinates = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(
                start=0, end=len(seq), seqname=source.subgraph_id,
                reading_frame_index=source.reading_frame_index
            )
        )
        var_node = self.create_node(
            seq=seq,
            variants=[variant_with_coordinates],
            reading_frame_index=source.reading_frame_index,
            orf=source.orf,
            level=self.root.level,
            global_variant=source.global_variant
        )

        if is_deletion:
            subgraph_id = self.subgraphs.generate_subgraph_id()
            var_node.subgraph_id = subgraph_id
            level = self.root.level + 1
            var_node.level = level
            self.subgraphs.add_subgraph(
                child_id=subgraph_id, parent_id=self.id, level=level,
                start=variant.location.start, end=variant.location.end,
                variant=variant, feature_type='transcript', feature_id=self.id
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

        # variant end
        if variant_end < target_end:
            index = target.seq.get_query_index(variant_end)
            head, tail = self.splice(target, index, 'reference')
            self.add_edge(var_node, tail, 'variant_end')
            returns[1] = head
            if in_frame:
                if variant_start == source_start:
                    returns[0] = head
                else:
                    returns[1] = returns[0]
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

    def insert_flanking_variant(self, cursors:List[TVGNode],
            var:VariantRecordWithCoordinate,
            seq:DNASeqRecordWithCoordinates,
            variants:List[VariantRecord], position:int,
            attach_directly:bool=False, known_orf_index:int=None,
            coordinate_feature_type:str=None, coordinate_feature_id:str=None
            ) -> List[TVGNode]:
        """ Insert flanking variant node into the graph. This function is used
        in apply_fusion to insert the intronic sequences.

        Args:
            cursors (List[TVGNode]): A list of three nodes that of the cursors
                during create_variant_graph.
            var (VariantRecordWithCoordinate): The variant object to be added
                to the inserted nodes. This should be the original variant to
                be applied, e.g. the fusion.
            seq (VariantRecordWithCoordinate): The sequence being inserted.
            variants (List[VariantRecord]): Additional variants that the
                insertion sequence carries.
            position (int): The position where the variant should be inserted.
                Only useful when `attach_directly` is false.
            attach_directly (bool): Whether attach the subgraph directly to
                cursors instead of checking for the position of insertion. This
                is useful to insert the fusion subgraph when an intronic
                subgraph is already inserted.
            known_orf_index (int): The known orf index (1, 2, or 3)
        """
        cursors = copy.copy(cursors)
        var_tails = []
        subgraph_id = self.subgraphs.generate_subgraph_id()
        for loc in seq.locations:
            loc.ref.seqname = subgraph_id
        branch = ThreeFrameTVG(
            seq, subgraph_id,
            global_variant=var.variant,
            cleavage_params=self.cleavage_params,
            max_adjacent_as_mnv=self.max_adjacent_as_mnv,
            coordinate_feature_type=coordinate_feature_type,
            coordinate_feature_id=coordinate_feature_id
        )
        level = cursors[0].level + 1
        branch.update_node_level(level)
        parent_id = cursors[0].subgraph_id
        if attach_directly:
            subgraph_start = cursors[0].seq.locations[-1].ref.end
            subgraph_end = subgraph_start + 1
        else:
            subgraph_start = var.variant.location.start
            subgraph_end = var.variant.location.end
        self.subgraphs.add_subgraph(
            child_id=subgraph_id, parent_id=parent_id, level=level,
            start=subgraph_start, end=subgraph_end,
            variant=var.variant,
            feature_type=coordinate_feature_type,
            feature_id=coordinate_feature_id
        )
        branch.init_three_frames(truncate_head=False)
        for rf_index, root in enumerate(branch.reading_frames):
            var_i = copy.deepcopy(var)
            var_i.location.reading_frame_index = rf_index
            var_i.location.seqname = subgraph_id
            node = list(root.out_edges)[0].out_node
            node.variants.append(var_i)
        branch.create_variant_graph(
            variants=variants,
            variant_pool=None,
            genome=None,
            anno=None,
            active_frames=[True, True, True],
            known_orf_index=known_orf_index,
            unmutated_start_size=0
        )

        for i in range(3):
            var_head = branch.reading_frames[i].get_reference_next()
            var_tail = var_head
            next_node = var_tail.get_reference_next()
            while next_node:
                var_tail = next_node
                next_node = var_tail.get_reference_next()
            var_tails.append(var_tail)
            while var_head.in_edges:
                edge = var_head.in_edges.pop()
                branch.remove_edge(edge)

            source = cursors[i]
            source_start = source.seq.locations[0].ref.start

            if attach_directly:
                self.add_edge(source, var_head, 'reference')
                continue

            if position < source_start:
                raise ValueError(
                    'The location of cursor is behind the variant. Something'
                    ' must have messed up during creating the variant graph.'
                )

            if position == source_start:
                prev = source.get_reference_prev()
                self.add_edge(prev, var_head, 'variant_start')
            else:
                index = source.seq.get_query_index(position)
                if index == -1:
                    if source.seq.locations[-1].ref.end != position:
                        raise ValueError(
                            'Fusion position is not in the range of the transcript'
                        )
                    self.add_edge(source, var_head, 'reference')
                else:
                    source, _ = self.splice(source, index, 'reference')
                    self.add_edge(source, var_head, 'variant_start')

        return var_tails

    def apply_fusion(self, cursors:List[TVGNode], variant:VariantRecord,
            variant_pool:VariantRecordPoolOnDisk, genome:DNASeqDict,
            anno:GenomicAnnotation,
            tx_seqs:Dict[str,DNASeqRecordWithCoordinates]=None,
            gene_seqs:Dict[str,DNASeqRecordWithCoordinates]=None,
            known_orf_index:int=None
            ) -> List[TVGNode]:
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
        original_cursors = cursors
        cursors = copy.copy(cursors)
        accepter_gene_id = variant.attrs['ACCEPTER_GENE_ID']
        accepter_tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
        accepter_tx_model = anno.transcripts[accepter_tx_id]
        accepter_chrom = accepter_tx_model.transcript.location.seqname
        exclude_variant_types = ['Fusion', 'Insertion', 'Deletion', 'Substitution', 'circRNA']
        if tx_seqs and accepter_tx_id in tx_seqs:
            accepter_tx_seq = tx_seqs[accepter_tx_id]
        else:
            accepter_tx_seq = accepter_tx_model.get_transcript_sequence(
                genome[accepter_chrom])

        # create subgraph for the left insertion
        if variant.attrs['LEFT_INSERTION_START'] is not None:
            tx_id = variant.location.seqname
            tx_model = anno.transcripts[tx_id]
            gene_id = tx_model.transcript.gene_id
            gene_model = anno.genes[gene_id]
            chrom = gene_model.chrom
            if gene_seqs and gene_id in gene_seqs:
                seq = gene_seqs[gene_id]
            else:
                seq = gene_model.get_gene_sequence(genome[chrom])
            insertion_start = variant.attrs['LEFT_INSERTION_START']
            insertion_end = variant.attrs['LEFT_INSERTION_END']
            insert_seq = seq[insertion_start:insertion_end]
            insertion_variants = variant_pool.filter_variants(
                gene_id=gene_id, exclude_type=exclude_variant_types,
                start=insertion_start, end=insertion_end, intron=True, return_coord='gene'
            )
            var = copy.deepcopy(variant)
            var.is_real_fusion = False
            var = seqvar.VariantRecordWithCoordinate(
                variant=var,
                location=FeatureLocation(start=0, end=len(insert_seq.seq))
            )
            cursors = self.insert_flanking_variant(
                cursors=cursors,
                var=var,
                seq=insert_seq,
                variants=insertion_variants,
                position=variant.location.start,
                attach_directly=False,
                known_orf_index=known_orf_index,
                coordinate_feature_type='gene',
                coordinate_feature_id=gene_id
            )

        # create subgraph for the right insertion
        if variant.attrs['RIGHT_INSERTION_START'] is not None:
            tx_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
            tx_model = anno.transcripts[tx_id]
            gene_id = tx_model.transcript.gene_id
            gene_model = anno.genes[gene_id]
            chrom = gene_model.chrom
            if gene_seqs and gene_id in gene_seqs:
                seq = gene_seqs[gene_id]
            else:
                seq = gene_model.get_gene_sequence(genome[chrom])
            insertion_start = variant.attrs['RIGHT_INSERTION_START']
            insertion_end = variant.attrs['RIGHT_INSERTION_END']
            insert_seq = seq[insertion_start:insertion_end]
            insertion_variants = variant_pool.filter_variants(
                gene_id=gene_id, exclude_type=exclude_variant_types,
                start=insertion_start, end=insertion_end, intron=True, return_coord='gene'
            )
            var = copy.deepcopy(variant)
            var.is_real_fusion = False
            var = seqvar.VariantRecordWithCoordinate(
                variant=var,
                location=FeatureLocation(start=0, end=len(insert_seq.seq))
            )
            if variant.attrs['LEFT_INSERTION_START'] is None:
                position = variant.location.start
                attach_directly = False
            else:
                position = None
                attach_directly = True
            cursors = self.insert_flanking_variant(
                cursors=cursors,
                var=var,
                seq=insert_seq,
                variants=insertion_variants,
                position=position,
                attach_directly=attach_directly,
                known_orf_index=known_orf_index,
                coordinate_feature_type='gene',
                coordinate_feature_id=gene_id
            )

        breakpoint_gene = variant.get_accepter_position()
        breakpoint_tx = anno.coordinate_gene_to_transcript(breakpoint_gene,
            accepter_gene_id, accepter_tx_id)
        accepter_seq = accepter_tx_seq[breakpoint_tx:]
        accepter_seq.orf = None
        accepter_variant_records = variant_pool.filter_variants(
            tx_ids=[accepter_tx_id], exclude_type=exclude_variant_types,
            start=breakpoint_gene, return_coord='transcript', intron=False
        )

        # add the variant as global variant
        accepter_length = len(accepter_tx_seq) - breakpoint_tx
        var = seqvar.VariantRecordWithCoordinate(
            variant=copy.deepcopy(variant),
            location=FeatureLocation(start=0, end=accepter_length)
        )
        var.variant.is_real_fusion = True

        if variant.attrs['LEFT_INSERTION_START'] is None and \
                variant.attrs['RIGHT_INSERTION_START'] is None:
            position = variant.location.start
            attach_directly = False
        else:
            position = None
            attach_directly = True

        self.insert_flanking_variant(
            cursors=cursors,
            var=var,
            seq=accepter_seq,
            variants=accepter_variant_records,
            position=position,
            attach_directly=attach_directly,
            known_orf_index=known_orf_index,
            coordinate_feature_type='transcript',
            coordinate_feature_id=accepter_tx_id
        )
        return original_cursors

    def _apply_insertion(self, cursors:List[TVGNode],
            var:seqvar.VariantRecordWithCoordinate,
            seq:DNASeqRecordWithCoordinates,
            variants:List[VariantRecord], active_frames:List[bool],
            coordinate_feature_type:str, coordinate_feature_id:str
        ) -> List[TVGNode]:
        """ This is a wrapper function to apply insertions to the graph. It can
        be used for actual Insertion, and also Substitution. It is also used
        in `apply_fusion` to to insert intronic sequences if the breakpoint is
        intronic.

        Args:
            cursors (List[TVGNode]): A list of three nodes that of the cursors
                during create_variant_graph.
            var (VariantRecordWithCoordinate): The variant object to be added
                to the inserted nodes.
            variants (List[VariantRecord]): Additional variants that the
                insertion sequence carries.
            active_frames (List[bool]): Whether each reading frame is active.
        """
        cursors = copy.copy(cursors)
        subgraph_id = self.subgraphs.generate_subgraph_id()
        for loc in seq.locations:
            loc.ref.seqname = subgraph_id
        branch = ThreeFrameTVG(
            seq, subgraph_id, global_variant=var.variant,
            max_adjacent_as_mnv=self.max_adjacent_as_mnv,
            coordinate_feature_type=coordinate_feature_type,
            coordinate_feature_id=coordinate_feature_id
        )
        level = cursors[0].level + 1
        branch.update_node_level(level)
        parent_id = cursors[0].subgraph_id
        self.subgraphs.add_subgraph(
            child_id=subgraph_id, parent_id=parent_id, level=level,
            start=var.location.start, end=var.location.end,
            variant=var.variant,
            feature_type=coordinate_feature_type,
            feature_id=coordinate_feature_id
        )
        branch.init_three_frames(truncate_head=False)
        for rf_index, root in enumerate(branch.reading_frames):
            var_i = copy.deepcopy(var)
            var_i.location.reading_frame_index = rf_index
            var_i.location.seqname = subgraph_id
            node = list(root.out_edges)[0].out_node
            node.variants.append(var_i)
            is_frameshifting = False
            if var.variant.is_frameshifting():
                is_frameshifting = True
        known_orf_index=self.get_known_reading_frame_index() \
                if self.has_known_orf else None
        branch.create_variant_graph(
            variants=variants,
            variant_pool=None,
            genome=None,
            anno=None,
            active_frames=active_frames,
            known_orf_index=known_orf_index
        )

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
                raise ValueError(
                    'The location of cursor is behind the variant. Something'
                    ' must have messed up during creating the variant graph.'
                )

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
            variant:VariantRecord,
            variant_pool:VariantRecordPoolOnDisk,
            genome:DNASeqDict, anno:GenomicAnnotation,
            gene_seqs:Dict[str, DNASeqRecordWithCoordinates]=None,
            active_frames:List[bool]=None
            ) -> List[TVGNode]:
        """ Apply an insertion into the the TVG graph. """
        gene_id = variant.attrs['DONOR_GENE_ID']
        gene_model = anno.genes[gene_id]
        chrom = gene_model.chrom
        donor_start = variant.get_donor_start()
        donor_end = variant.get_donor_end()
        if gene_seqs and gene_id in gene_seqs:
            gene_seq = gene_seqs[gene_id]
        else:
            gene_seq = gene_model.get_gene_sequence(genome[chrom])
        insert_seq = gene_seq[donor_start:donor_end]
        exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion', 'circRNA']
        insert_variants = variant_pool.filter_variants(
            gene_id=gene_id, exclude_type=exclude_type,
            start=donor_start, end=donor_end
        )
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
        if variant.is_end_inclusion():
            insert_seq = insert_seq + ref_seq
        else:
            insert_seq = ref_seq + insert_seq
        var = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=1, end=len(insert_seq.seq))
        )
        return self._apply_insertion(
            cursors=cursors,
            var=var,
            seq=insert_seq,
            variants=insert_variants,
            active_frames=active_frames,
            coordinate_feature_type='gene',
            coordinate_feature_id=gene_id
        )

    def apply_substitution(self, cursors:List[TVGNode],
            variant:VariantRecord,
            variant_pool:VariantRecordPoolOnDisk,
            genome:DNASeqDict, anno:GenomicAnnotation,
            gene_seqs:Dict[str, DNASeqRecordWithCoordinates]=None,
            active_frames:List[bool]=None
            ) -> List[TVGNode]:
        """ Apply a substitution variant into the graph """
        gene_id = variant.attrs['DONOR_GENE_ID']
        gene_model = anno.genes[gene_id]
        chrom = gene_model.chrom
        donor_start = variant.get_donor_start()
        donor_end = variant.get_donor_end()
        if gene_seqs and gene_id in gene_seqs:
            gene_seq = gene_seqs[gene_id]
        else:
            gene_seq = gene_model.get_gene_sequence(genome[chrom])
        sub_seq = gene_seq[donor_start:donor_end]
        exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion', 'circRNA']
        sub_variants = variant_pool.filter_variants(
            gene_id=gene_id, exclude_type=exclude_type,
            start=donor_start, end=donor_end
        )
        var = seqvar.VariantRecordWithCoordinate(
            variant=variant,
            location=FeatureLocation(start=0, end=len(sub_seq.seq))
        )
        return self._apply_insertion(
            cursors=cursors,
            var=var,
            seq=sub_seq,
            variants=sub_variants,
            active_frames=active_frames,
            coordinate_feature_type='gene',
            coordinate_feature_id=gene_id
        )


    def create_variant_graph(self, variants:List[VariantRecord],
            variant_pool:Union[VariantRecordWithCoordinate, VariantRecordPoolOnDisk],
            genome:DNASeqDict, anno:GenomicAnnotation,
            tx_seqs:Dict[str, DNASeqRecordWithCoordinates]=None,
            gene_seqs:Dict[str, DNASeqRecordWithCoordinates]=None,
            active_frames:List[bool]=None, known_orf_index:int=None,
            unmutated_start_size:int=3) -> None:
        """ Create a variant graph.

        With a list of genomic variants, incorprate each variant into the
        graph.

        For protein coding transcripts, in frame variants are not incorporated
        into non-canonical reading frames until the first frameshifting mutate
        shows up.

        Args:
            variant [seqvar.VariantRecord]: The variant record.
            unmutated_start_size [int]: Number of nucleotides after the cds start
                site (or 0 for noncoding) that should not contain any variant.
        """
        ## Fitler variants
        # skipping start lost mutations
        start_index = self.seq.orf.start if self.has_known_orf else 0
        start_index += unmutated_start_size
        is_fusion = any(v.is_fusion() for v in variants)

        filtered_variants = []
        for v in variants:
            if v.location.start == start_index - 1 \
                    and (v.is_insertion() or v.is_deletion()) \
                    and not v.is_fusion() \
                    and not v.is_alternative_splicing():
                v.to_end_inclusion(self.seq)

            # Skip variants that the position is smaller than the first NT
            # after start codon. Exception for fusion, that if the donor
            # breakpoint is at the last NT of the start codon it is retained
            # because it won't break the start codon.
            if v.location.start < start_index and not \
                    (v.is_fusion() and v.location.start == start_index - 1):
                continue

            # if the transcript is mrna_end_NF, we are not going to use any
            # variants in the annotated 3'UTR region.
            if self.mrna_end_nf and not is_fusion:
                tx_end = self.seq.orf.end if self.seq.orf is not None \
                    else len(self.seq.seq)
                orf_end_trinuc = FeatureLocation(
                    start=tx_end - 3, end=tx_end
                )
                if v.location.overlaps(orf_end_trinuc):
                    continue

            filtered_variants.append(v)

        merged_mnvs = seqvar.find_mnvs_from_adjacent_variants(
            filtered_variants, self.max_adjacent_as_mnv
        )
        variants_with_mnv = sorted(filtered_variants + merged_mnvs)
        variant_iter = iter(variants_with_mnv)
        variant = next(variant_iter, None)
        cursors = copy.copy([x.get_reference_next() for x in self.reading_frames])

        if active_frames is None:
            if self.has_known_orf:
                if known_orf_index is None:
                    known_orf_index = self.get_known_reading_frame_index()
                active_frames = [False, False, False]
                active_frames[known_orf_index] = True
            else:
                active_frames = [True, True, True]

        while variant:
            if not any(cursors):
                break

            if any(c.seq.locations[0].ref.start > variant.location.start for c in cursors):
                variant = next(variant_iter, None)
                continue

            if variant.is_frameshifting():
                for i in range(3):
                    if not active_frames[i]:
                        continue
                    frames_shifted = variant.frames_shifted()
                    j = (i + frames_shifted) % 3
                    active_frames[j] = True

            # if any(not x for x in active_frames) and variant.is_frameshifting():
            #     frames_shifted = variant.frames_shifted()
            #     i = (known_orf_index + frames_shifted) % 3
            #     active_frames[i] = True

            any_cursor_expired = False
            for i,cursor in enumerate(cursors):
                if cursor.seq.locations[-1].ref.end == variant.location.start and \
                        variant.is_fusion() and variant.location.start == len(self.seq):
                    continue
                if cursor.seq.locations[-1].ref.end <= variant.location.start:
                    cursors[i] = cursors[i].get_reference_next()
                    any_cursor_expired = True
            if any_cursor_expired:
                continue

            if variant.type == 'Fusion':
                # for fusion, whether there is any frameshifting mutation does
                # not affect the main sequence, so a copy of active frames
                # are parsed to the apply_fusion function.
                cursors = self.apply_fusion(
                    cursors=cursors,
                    variant=variant,
                    variant_pool=variant_pool,
                    genome=genome,
                    anno=anno,
                    tx_seqs=tx_seqs,
                    gene_seqs=gene_seqs,
                    known_orf_index=known_orf_index
                )

            elif variant.type == 'Insertion':
                cursors = self.apply_insertion(
                    cursors=cursors,
                    variant=variant,
                    variant_pool=variant_pool,
                    genome=genome,
                    anno=anno,
                    gene_seqs=gene_seqs,
                    active_frames=active_frames
                )

            elif variant.type == 'Substitution':
                cursors = self.apply_substitution(
                    cursors=cursors,
                    variant=variant,
                    variant_pool=variant_pool,
                    genome=genome,
                    anno=anno,
                    gene_seqs=gene_seqs,
                    active_frames=active_frames
                )

            elif variant.is_frameshifting():
                frames_shifted = variant.frames_shifted()
                for i in range(3):
                    if not active_frames[i]:
                        continue
                    j = (i + frames_shifted) % 3
                    cursors[i], cursors[j] = self.apply_variant(cursors[i],
                        cursors[j], variant)
            else:
                for i,_ in enumerate(cursors):
                    if not active_frames[i]:
                        continue
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
        node_copy = node.copy()
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

    def merge_into_inbonds(self, node:TVGNode) -> List[TVGNode]:
        """ Merge the target node into its inbond nodes. """
        in_nodes = node.get_in_nodes()
        for in_node in in_nodes:
            if len(in_node.out_edges) > 1:
                new_in_node = in_node.copy()
                for edge in in_node.in_edges:
                    upstream = edge.in_node
                    self.add_edge(upstream, new_in_node, edge.type)
                new_in_node.append_right(node)
                for edge in node.out_edges:
                    self.add_edge(new_in_node, edge.out_node, edge.type)
            else:
                in_node.append_right(node)
                in_node.out_edges = set()
                for edge in node.out_edges:
                    self.add_edge(in_node, edge.out_node, edge.type)
        self.remove_node(node)
        return in_nodes

    def find_bridge_nodes_between(self, start:TVGNode, end:TVGNode, members:Set[TVGNode]
            ) -> Tuple[Set[TVGNode], Set[TVGNode], Set[TVGNode], Set[TVGNode]]:
        """ Find all bridge in and bridge out nodes between a start and a
        end node. """
        this_id = start.get_last_rf_index()
        bridge_in:Set[TVGNode] = set()
        bridge_out:Set[TVGNode] = set()
        subgraph_in:Set[TVGNode] = set()
        subgraph_out:Set[TVGNode] = set()
        visited = set()
        queue:Deque[TVGNode] = deque([e.out_node for e in start.out_edges])

        # This `end_nodes` is to set the boundray of the current variant bubble.
        # Any "inbrige" and "outbrange" nodes to and from any node outside of
        # the current variant bubble should not be considered.
        end_nodes = {}
        for n in members:
            for out_node in n.get_out_nodes():
                if out_node not in members:
                    end_nodes[out_node.id] = out_node

        while queue:
            cur = queue.pop()
            len_before = len(visited)
            visited.add(cur)
            len_after = len(visited)
            if len_before == len_after:
                continue

            is_bridge_out = any(
                e.out_node.get_first_rf_index() != this_id
                    or e.out_node.get_last_rf_index() != this_id
                    and any(
                        v.variant.is_frameshifting()
                        for v in cur.variants
                    )
                    and not cur.is_reference()
                    and e.out_node not in members
                for e in cur.out_edges
            )
            # hybrid bridge out is when the node is already processed at least
            # once by merging it with downstream nodes, so different segment.
            # of the node has different rf index. This is to handle a combination
            # of indel on AltSplice Ins/Sub that goes back to original reading
            # frame. We may have to change it to check whether any segment of
            # the node has different rf index. Issue #726
            is_hybrid_bridge_out = cur.get_last_rf_index() != this_id \
                and not cur.is_reference()

            is_internal_bridge = all(x in members for x in cur.get_in_nodes()) \
                and all(x in members for x in cur.get_out_nodes())

            is_bridge_out |= (is_hybrid_bridge_out and not is_internal_bridge)

            if is_bridge_out and cur is not end:
                bridge_out.add(cur)
                continue

            if not self.is_circ_rna():
                if cur not in members:
                    if cur.subgraph_id != start.subgraph_id:
                        if not cur.is_inframe_subgraph(start, end):
                            subgraph_out.add(cur)
                            continue

                    if cur.has_multiple_segments() \
                            and cur.get_first_subgraph_id() == start.subgraph_id \
                            and cur.get_first_subgraph_id() != cur.get_last_subgraph_id():
                        subgraph_out.add(cur)
                        continue
                    if cur.id in end_nodes:
                        continue
                elif cur is not end:
                    # This is when an altSplice indel/sub carries additional frameshift
                    # variant, and the indel frameshift is very closed to the end of
                    # the subgraph, there will be additional subgraph end nodes
                    # that go back to the main graph.
                    max_level = self.subgraphs[cur.get_max_subgraph_id(self.subgraphs)].level
                    if all(max_level > self.subgraphs[n.subgraph_id].level and n not in members
                            for n in cur.get_out_nodes()):
                        subgraph_out.add(cur)

            for e in cur.in_edges:
                if e.in_node.get_first_rf_index() != this_id or e.in_node.was_bridge \
                        and e.in_node not in visited and e.in_node is not start:
                    bridge_in.add(e.in_node)
                elif not self.is_circ_rna() \
                        and e.in_node.subgraph_id != cur.subgraph_id:
                    if e.in_node.subgraph_id != start.subgraph_id:
                        subgraph_in.add(e.in_node)
                    elif e.in_node.is_inframe_subgraph(start, end):
                        if not e.in_node in members:
                            subgraph_in.add(e.in_node)

            if cur is not end:
                for e in cur.out_edges:
                    queue.appendleft(e.out_node)
        for node in copy.copy(subgraph_in):
            if node in subgraph_out:
                subgraph_in.remove(node)
                subgraph_out.remove(node)

        inframe_singleton_subgraphs = set()
        for node in subgraph_out:
            if node in subgraph_in:
                inframe_singleton_subgraphs.add(node)

        subgraph_in = {x for x in subgraph_in if x not in inframe_singleton_subgraphs}
        subgraph_out = {x for x in subgraph_out if x not in inframe_singleton_subgraphs}

        return bridge_in, bridge_out, subgraph_in, subgraph_out

    @staticmethod
    def find_other_subgraph_end_nodes(node:TVGNode, members:Set[TVGNode]) -> Set[TVGNode]:
        """ Find parallel subgraph end nodes """
        end_nodes = set()
        for out_node in node.get_out_nodes():
            if out_node.level >= node.level:
                raise ValueError('something went wrong!')
            for end_node in out_node.get_in_nodes():
                if end_node is not node and \
                        end_node.subgraph_id == node.subgraph_id and\
                        end_node.reading_frame_index == node.reading_frame_index:
                    end_nodes.add(end_node)
        for in_node in node.get_in_nodes():
            for out_node in in_node.get_out_nodes():
                if not any(x in members for x in out_node.get_out_nodes()) \
                        and out_node.get_first_rf_index() == node.get_first_rf_index() \
                        and out_node.get_first_subgraph_id() == node.get_first_subgraph_id():
                    end_nodes.add(out_node)
        return end_nodes

    @staticmethod
    def find_other_member_end_nodes(members:Iterable[TVGNode]) -> Set[TVGNode]:
        """ Find other member end nodes """
        end_nodes = set()
        for member in members:
            for out_node in member.get_out_nodes():
                if out_node not in members:
                    end_nodes.add(member)
        return end_nodes

    def is_fusion_subgraph_out(self, up:TVGNode, down:TVGNode) -> bool:
        """ Checks if the given node is a fusion subgraph out """
        if not up.is_inbond_of(down):
            return False

        if not any(v.variant.is_fusion() for v in down.variants):
            return False

        up_subgraph = self.subgraphs[up.subgraph_id]
        down_subgraph = self.subgraphs[down.subgraph_id]

        return up_subgraph.is_parent(down_subgraph)

    def find_variant_bubble(self, node:TVGNode, min_size:int=6
            ) -> Tuple[TVGNode, Set[TVGNode]]:
        r""" Find the variable bubble, and return the end node and all
        member nodes of the bubble. The end node is defined as the first node
        downstream that no any variant is at its location, and the length is
        at longer or equal to `min_size`.

        For example, in a graph like below, the node ATGG's farthest node
        with overlap would be node 'CCCT'
                 T--
                /   \
            ATGG-TCT-G-CCCT
                    \ /
                     A
        """
        # When the subgraph_checker is True, alignment is limited within the
        # current subgraph. The only exception is for fusions, that incoming
        # node in the main graph is exclusively bond with the fusion donor.
        subgraph_checker = True
        if self.is_circ_rna():
            subgraph_checker = False
        # find the range of overlaps
        farthest = None
        if not node.get_reference_next():
            return None, set()
        if not node.get_reference_next().out_edges:
            return node.get_reference_next(), set()

        start_variants = {v.variant for v in node.variants}

        def is_candidate_out_node(x:TVGNode, y:TVGNode):
            # Note: have to use y.subgraph_id because for deletion, subgraph_id
            # is different from get_first_subgraph_id()
            # try:
            #     is_in_subgraph = x.get_last_subgraph_id() == y.subgraph_id
            # except IndexError:
            #     is_in_subgraph = x.subgraph_id == y.subgraph_id
            return subgraph_checker is False \
                or x.get_max_subgraph_id(self.subgraphs) == y.subgraph_id \
                or self.is_fusion_subgraph_out(x,y)

        queue:Deque[TVGNode] = deque([])
        for out_node in node.get_out_nodes():
            if is_candidate_out_node(node, out_node):
                queue.append(out_node)

        if not queue:
            # This means the input node is an end of a subgraph of alt splice
            # insertion/deletion, so no further aligning is needed. Thus return
            # the input node itself.
            return node, set()

        visited:Dict[str, TVGNode] = {node.id: node}
        # When a new farthest node is found, the downstream nodes of the old
        # farthest node are put into this `exceptions` container, so they
        # can be visited again.
        exceptions:Set[TVGNode] = set()
        non_members:Set[str] = set()
        while queue:
            cur:TVGNode = queue.popleft()
            if cur is None:
                continue

            if cur.reading_frame_index != node.reading_frame_index:
                non_members.add(cur.id)
                continue

            cur_variants = {v.variant for v in cur.variants}
            if any(v.is_frameshifting() and v not in start_variants for v in cur_variants):
                continue

            if subgraph_checker:
                if cur.has_exclusive_inbond_node() and \
                        cur.get_in_nodes()[0].subgraph_id == node.subgraph_id:
                    subgraph_checker = False

            visited_len_before = len(visited)
            visited[cur.id] = cur
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                if cur is farthest and cur is not node:
                    if not cur.out_edges:
                        continue

                    downstream = cur.get_reference_next()
                    downstream_has_subgraph_in = \
                        not self.is_circ_rna() \
                        and not any(x.variant.is_fusion() for x in downstream.variants) \
                        and downstream.has_in_subgraph()

                    if cur.get_reference_next().has_in_bridge() \
                            or downstream_has_subgraph_in:
                        for out_node in cur.get_out_nodes():
                            if is_candidate_out_node(cur, out_node):
                                queue.append(out_node)
                    # if the farthest has less than 6 neucleotides, continue
                    # searching, because it's likely not able to keep one amino
                    # acid after translation.
                    elif len(cur.seq) < min_size:
                        for out_node in cur.get_out_nodes():
                            if is_candidate_out_node(cur, out_node):
                                queue.append(out_node)
                elif cur in exceptions:
                    for out_node in cur.get_out_nodes():
                        if is_candidate_out_node(cur, out_node):
                            queue.append(out_node)
                            exceptions.add(out_node)
                    exceptions.remove(cur)
                continue

            if farthest is None and cur.is_reference():
                farthest = cur
                queue.append(cur)
                continue

            if not cur.is_reference() \
                    or ( len(cur.out_edges) > 1
                        and any(len(x.in_edges) > 1 for x in cur.get_out_nodes())):
                for out_node in cur.get_out_nodes():
                    if is_candidate_out_node(cur, out_node):
                        queue.append(out_node)
                continue

            if cur.is_stop_node():
                continue

            if self.first_node_is_smaller(cur, farthest):
                for out_node in cur.get_out_nodes():
                    if is_candidate_out_node(cur, out_node):
                        queue.append(out_node)
                continue

            farthest, cur = cur, farthest
            queue.append(cur)
            queue.append(farthest)
            exceptions.add(cur)
            continue
        members = {v for k,v in visited.items() if k not in non_members}
        return farthest, members

    def first_node_is_smaller(self, first:TVGNode, second:TVGNode) -> bool:
        """ Check if the first node is larger """
        if first.subgraph_id == second.subgraph_id:
            return first.seq.locations[0] < second.seq.locations[0]

        subgraph1, subgraph2 = self.subgraphs.find_compatible_parents(first, second)

        if subgraph1.id == subgraph2.id:
            raise ValueError('They should not equal. Something went wrong.')

        if subgraph1.id == subgraph2.parent_id:
            if subgraph1.id != first.subgraph_id:
                raise ValueError('They should not equal. Something went wrong.')
            return first.seq.locations[0].ref.start < subgraph2.location.end

        if subgraph2.id == subgraph1.parent_id:
            if subgraph2.id != second.subgraph_id:
                raise ValueError('They should not equal. Something went wrong.')
            return subgraph1.location.end < second.seq.locations[0].ref.start

        return subgraph1.location < subgraph2.location

    @staticmethod
    def nodes_have_too_many_variants(nodes:Iterable[TVGNode],
            max_in_bubble_variants:int) -> bool:
        """ Check the total number of variants of given nodes """
        if max_in_bubble_variants == -1:
            return False
        variants = set()
        for node in nodes:
            for variant in node.variants:
                variants.add(variant.variant)
        return len(variants) > max_in_bubble_variants

    @staticmethod
    def get_max_in_bubble_variants(n:int) -> int:
        """ Get the `max_in_bubble_variants` based on the total number of variants
        in a variant bubble. The values are set so the number of combinations to
        consider won't be too much more than 5,000. """
        if n <= 12:
            return -1
        if n <= 14:
            return 7
        if n <= 15:
            return 6
        if n <= 16:
            return 5
        if n <= 21:
            return 4
        if n <= 34:
            return 3
        if n <= 100:
            return 2
        return 1

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
                be aligned
        Returns:
            The original input node.
        """
        start_node = node
        rf_index = node.reading_frame_index
        if rf_index is None:
            raise ValueError('reading_frame_index not found')
        end_node, members = self.find_variant_bubble(node)

        member_variants = set()
        for member in members:
            member_variants.update([v.variant.id for v in member.variants])

        max_in_bubble_variants = self.get_max_in_bubble_variants(len(member_variants))

        if len(member_variants) >= 13:
            get_logger().warning(
                "Hypermutated region detected with %i variants."
                " `max_in_bubble_variants` of %i is used.",
                len(member_variants), max_in_bubble_variants
            )

        if not end_node:
            raise err.FailedToFindVariantBubbleError()

        end_nodes = {end_node.id}
        if not self.is_circ_rna() and end_node.is_subgraph_end():
            subgraph_ends = {end_node}
            subgraph_ends.update(self.find_other_subgraph_end_nodes(end_node, members))
            if len(subgraph_ends) > 1:
                for subgraph_end in subgraph_ends:
                    end_nodes.update([x.id for x in subgraph_end.get_out_nodes()])

        bridges = self.find_bridge_nodes_between(start_node, end_node, members)
        bridge_ins, bridge_outs, subgraph_ins, subgraph_outs = bridges

        for bridge in bridge_outs:
            for e in bridge.out_edges:
                end_nodes.add(e.out_node.id)

        for subgraph in subgraph_outs:
            for e in subgraph.out_edges:
                if e.out_node not in members:
                    end_nodes.add(e.out_node.id)

        new_nodes:Set[TVGNode] = set()
        queue = deque()
        trash = set()
        bridge_map:Dict[TVGNode, TVGNode] = {}
        for bridge_in in list(bridge_ins) + list(subgraph_ins):
            new_bridge = bridge_in.copy()
            for edge in bridge_in.out_edges:
                self.add_edge(new_bridge, edge.out_node, edge.type)
                # Here we want to limit the process within the variant bubble.
                # In-bridge nodes should not be merged with their outgoing nodes
                # that do not belong to the bubble.
                if edge.out_node not in members:
                    end_nodes.add(edge.out_node.id)
            bridge_map[new_bridge] = bridge_in
            trash.add(bridge_in)

        new_bridges:Set[TVGNode] = set()

        # start by removing the out edges from the start node.
        for edge in copy.copy(start_node.out_edges):
            out_node = edge.out_node
            if out_node in subgraph_outs:
                continue
            self.remove_edge(edge)
            queue.appendleft(out_node)

        for bridge in bridge_map:
            queue.appendleft(bridge)

        while queue:
            cur:TVGNode = queue.pop()
            if cur.id in end_nodes or not cur.out_edges:
                if cur not in bridge_map:
                    new_nodes.add(cur)
                else:
                    new_bridges.add(cur)
                continue

            # When all out nodes are in `end_nodes`. This is to avoid the `cur`
            # to be replicated.
            if all(x.id in end_nodes for x in cur.get_out_nodes()):
                new_node = cur.copy()
                trash.add(cur)
                for edge in cur.out_edges:
                    self.add_edge(new_node, edge.out_node, _type=edge.type)
                if cur not in bridge_map:
                    new_nodes.add(new_node)
                else:
                    bridge_map[new_node] = bridge_map[cur]
                    new_bridges.add(new_node)
                continue

            for out_edge in copy.copy(cur.out_edges):
                out_node:TVGNode = out_edge.out_node

                # So this is the case that some of the out_nodes are in end_nodes
                # but not the others.
                if out_node.id in end_nodes:
                    new_node = cur.copy()
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

                if self.nodes_have_too_many_variants([cur, out_node], max_in_bubble_variants):
                    continue

                # create new node with the combined sequence
                new_node = cur.copy()
                was_bridge = new_node.reading_frame_index != rf_index \
                    and any(x.reading_frame_index != rf_index for x in out_node.get_out_nodes())
                new_node.append_right(out_node)
                new_node.was_bridge = was_bridge
                if new_node.level < out_node.level:
                    new_node.subgraph_id = out_node.subgraph_id
                    new_node.level = out_node.level

                for edge in copy.copy(out_node.out_edges):
                    edge_type = 'reference' \
                        if self.is_reference_edge(new_node, edge.out_node) \
                        else 'variant_end'
                    self.add_edge(new_node, edge.out_node, _type=edge_type)

                if out_node.id not in end_nodes:
                    queue.appendleft(new_node)

                if cur in bridge_map:
                    bridge_map[new_node] = bridge_map[cur]
            # now remove the cur node from graph
            if cur not in bridge_map:
                self.remove_node(cur)

        # add new nodes to the graph
        for new_node in new_nodes:
            if self.is_reference_edge(start_node, new_node):
                in_edge_type = 'reference'
            else:
                in_edge_type = 'variant_start'
            self.add_edge(start_node, new_node, in_edge_type)

        for bridge in new_bridges:
            original_bridge = bridge_map[bridge]
            for edge in original_bridge.in_edges:
                self.add_edge(edge.in_node, bridge, edge.type)
            trash.add(original_bridge)

        for trash_node in trash:
            self.remove_node(trash_node)

    def expand_alignments(self, start:TVGNode) -> List[TVGNode]:
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
        downstream_empty_error_msg = 'Downstream node becomes empty when it is' \
            + 'not the end of the graph. Something is wrong.'
        if len(start.out_edges) == 1 and len(start.get_out_nodes()[0].in_edges) > 1:
            downstream = start.get_out_nodes()[0]
            right_index = (3 - len(start.seq) % 3) % 3
            if right_index >= len(downstream.seq.seq):
                self.merge_into_inbonds(downstream)
                # The downstream should be end of a subgraph so nothing needs
                # to be returned.
                return []
            right_over = downstream.truncate_left(right_index)
            for in_node in downstream.get_in_nodes():
                in_node.append_right(right_over)
            if downstream.seq.seq == '':
                if downstream.get_out_nodes():
                    raise ValueError(downstream_empty_error_msg)
                self.remove_node(downstream)
                return []
            if downstream.is_subgraph_end():
                return []
            return [downstream]
        # number of NT to be carried over to the downstreams.
        if start.seq:
            left_index = len(start.seq) - len(start.seq) % 3
            left_over = start.truncate_right(left_index)
        else:
            left_over_seq = dna.DNASeqRecordWithCoordinates(Seq(''), [])
            left_over = self.create_node(
                seq=left_over_seq,
                subgraph_id=start.subgraph_id,
                level=start.level
            )

        end_nodes:List[TVGNode] = []

        for out_edge in start.out_edges:
            out_node = out_edge.out_node
            out_node.append_left(left_over)

        ref_node = start.get_reference_next()
        for edge in start.out_edges:
            out_node = edge.out_node
            if start.subgraph_id != out_node.subgraph_id and\
                    not any(x.out_node.subgraph_id == start.subgraph_id \
                        for x in out_node.out_edges) \
                    and not out_node is ref_node \
                    and not out_node.is_bridge() \
                    and not set(out_node.get_out_nodes()) == set(ref_node.get_out_nodes()):
                end_nodes.append(out_node)

        if not ref_node.out_edges:
            return end_nodes

        if len(ref_node.out_edges) > 1:
            end_nodes.append(ref_node)
            return end_nodes

        end = ref_node.get_reference_next()
        right_index = (3 - len(ref_node.seq) % 3) % 3

        if right_index >= len(end.seq.seq):
            if len(end.get_out_nodes()) == 1 \
                    and len(end.get_out_nodes()[0].get_in_nodes()) == 1:
                # This is when the left or right intronic insertion of a fusion
                # is smaller than 3. The `end` should contain an unique out node
                end = self.merge_with_outbonds(end)[0]
            elif end.global_variant and end.global_variant.is_fusion() \
                    and ref_node.has_exclusive_outbond_node():
                end = self.merge_with_outbonds(ref_node)[0]
                return [end]
            else:
                # Similar to above but here for AltSplice.
                self.merge_into_inbonds(end)
                return []

        right_over = end.truncate_left(right_index)

        for in_edge in end.in_edges:
            in_node = in_edge.in_node
            in_node.append_right(right_over)

        if end.seq.seq == '':
            if end.get_out_nodes():
                raise ValueError(downstream_empty_error_msg)
            self.remove_node(end)

        if end.subgraph_id != self.id and any(x.out_node.subgraph_id == self.id
                for x in end.out_edges) \
                or ref_node.subgraph_id != self.id and end.subgraph_id == self.id:
            return end_nodes
        end_nodes.append(end)
        return end_nodes

    def collapse_equivalent_nodes(self, node:TVGNode):
        """ Collapse nodes that have the same sequence and same incoming and
        outgoing nodes. """
        group:Dict[Tuple(Tuple[PVGNode], Tuple[PVGNode]),TVGNodeCollapser] = {}
        nodes = node.get_out_nodes()
        for out_node in nodes:
            in_nodes = tuple(out_node.get_in_nodes())
            out_nodes = tuple(out_node.get_out_nodes())
            collapser = group.setdefault((in_nodes, out_nodes), TVGNodeCollapser())
            redundant_node = collapser.collapse(out_node)
            if redundant_node:
                self.remove_node(redundant_node)
        collapsed_nodes = set(nodes)
        for _node in copy.copy(collapsed_nodes):
            if _node.is_orphan():
                collapsed_nodes.remove(_node)
        return collapsed_nodes

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

            try:
                self.align_variants(cur)
            except err.FailedToFindVariantBubbleError:
                continue

            self.collapse_equivalent_nodes(cur)
            if cur.out_edges:
                end_nodes = self.expand_alignments(cur)
                for node in end_nodes:
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
        root = PVGNode(None, None, subgraph_id=self.id, level=self.root.level)
        if self.has_known_orf:
            known_orf = [int(self.seq.orf.start), int(self.seq.orf.end)]
        else:
            known_orf = [None, None]

        pgraph = PeptideVariantGraph(
            root=root,
            _id=self.id,
            known_orf=known_orf,
            cleavage_params=self.cleavage_params,
            cds_start_nf=self.cds_start_nf,
            hypermutated_region_warned=self.hypermutated_region_warned,
            global_variant=self.global_variant,
            subgraphs=self.subgraphs,
            gene_id=self.gene_id,
        )

        queue = deque([(dnode, root) for dnode in self.reading_frames])
        visited = {}
        terminal_nodes:List[Tuple[PVGNode, int]] = []

        while queue:
            dnode, pnode = queue.pop()
            if not dnode.out_edges:
                pgraph.add_stop(pnode)
                if self.should_clip_trailing_nodes():
                    pnode.truncated = True
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

                out_node.check_stop_altering(self.seq.seq, orf[1])

                new_pnode = out_node.translate()

                if not self.is_circ_rna():
                    new_pnode.fix_selenocysteines(self.sect_variants, self.subgraphs)

                if new_pnode.seq.seq == '' and not out_node.get_out_nodes():
                    if not self.should_clip_trailing_nodes():
                        new_pnode.seq.seq = Seq('*')
                    elif len(dnode.out_edges) == 1:
                        pnode.truncated = True

                if orf[1] and out_node.level == 0:
                    orf_end_query = out_node.seq.get_query_index(orf[1])
                    if 0 < orf_end_query <= len(out_node.seq.seq) - 3 \
                            and orf_end_query % 3 == 0:
                        pnode_orf_end = int(orf_end_query / 3)
                        is_stop_lost = any(orf_end_query in v.location for v in new_pnode.variants)
                        if new_pnode.seq.seq[pnode_orf_end] != '*' \
                                and not is_stop_lost:
                            terminal_nodes.append((new_pnode, pnode_orf_end))

                new_pnode.orf = orf
                pnode.add_out_edge(new_pnode)
                visited[out_node] = new_pnode
                #if not hit_stop:
                queue.appendleft((out_node, new_pnode))

        for i, dnode in enumerate(self.reading_frames):
            dnode = dnode.get_reference_next()
            if dnode:
                pgraph.reading_frames[i] = visited[dnode]

        # Sometimes the annotated cds end isn't in fact a stop codon. If so,
        # a fake stop codon is added so the node won't be called as variant
        # peptide.
        if self.has_known_orf:
            for terminal_node, stop_site in terminal_nodes:
                right_part = terminal_node.split_node(stop_site)
                if(len(right_part.seq.seq) > 1 and right_part.seq.seq[0] == '*'):
                    right_part.split_node(1)
                else:
                    fake_stop = AminoAcidSeqRecordWithCoordinates(
                        seq='*', locations=[]
                    )
                    fake_stop_node = PVGNode(
                        seq=fake_stop,
                        reading_frame_index=terminal_node.reading_frame_index,
                        subgraph_id=self.id
                    )
                    terminal_node.add_out_edge(fake_stop_node)

        return pgraph

    def jsonfy(self):
        """ Create node and edge list from a ThreeFrameTVG object. """
        queue = deque([self.root])
        node_index:dict[int, int] = {}
        cur_index = 0

        nodes = []
        while queue:
            cur = queue.pop()
            if id(cur) in node_index:
                continue

            node = {
                'index': cur_index
            }
            node_index[id(cur)] = cur_index


            node['seq'] = str(cur.seq.seq) if cur.seq  else ''
            variants = set()
            for v in cur.variants:
                if v.variant.attrs.get('MERGED_MNV'):
                    variants.update(v.variant.attrs['INDIVIDUAL_VARIANT_IDS'])
                else:
                    variants.add(v.variant.id)
            node['variants'] = list(variants)
            node['rf_index'] = cur.reading_frame_index

            if self.has_known_orf \
                    and cur.seq \
                    and self.get_known_reading_frame_index() \
                        == cur.reading_frame_index:
                node['start_codon'] = cur.seq.get_query_index(
                    ref_index=self.seq.orf.start
                )
            else:
                node['start_codon'] = -1

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
