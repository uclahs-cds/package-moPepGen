""" Module for node in the peptide graph """
from __future__ import annotations
import copy
from functools import cmp_to_key
import math
from typing import Dict, List, Set, Iterable
from moPepGen import aa, circ, seqvar, err
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.svgraph.SubgraphTree import SubgraphTree


class PVGNode():
    """ The PVGNode class is used in the svgraph.PeptideVariantGraph, that
    stores the sequence in a Node.

    Attributes:
        - `seq` (AminoAcidSeqRecord): The amino acid sequence.
        - `reading_frame_index` (int): Reading frame index that the node belongs
          to. Must be 0, 1, or 2.
        - `variants` (List[seqvar.VariantRecordWithCoordinate]): The variant
          records carried in the sequence.
        - `in_nodes` (Set[PVGNode]): Inbound nodes.
        - `out_nodes` (Set[PVGNode]): Outbound nodes
        - `cleavage` (bool): Whether the start of the node is a cleavage site.
        - `truncated` (bool): Whether the node is truncated. Useful when the
          sequence does not have a confirmed stop codon.
        - `orf` (List[int]): List of two integers indicating the start and stop
          of the orf.
        - `was_bridge` (bool): When true indicating that the node used to be a
          bridge node between reading frames, but is no longer a bridge node
          because it was split.
        - `pre_cleaved` (bool): Indicating that the node is already cleaved.
        - `level` (bool): Subgraph level that the node belongs to.
        - `npop_collapsed` (bool): Whether the n-terminus is pop collapsed.
        - `cpop_collapsed` (bool): Whether the c-terminus is pop collapsed.
        - `upstream_indel_map` (Dict[PVGNode, List[VariantRecord]]): A mapping
          from upstream node to indel variants. When a node contains INDELs is
          collapsed with a canonical peptide node, the upstream node will be
          called for miscleaved peptides that has corresponding INDELs.
    """
    def __init__(self, seq:aa.AminoAcidSeqRecordWithCoordinates,
            reading_frame_index:int, subgraph_id:str,
            variants:List[seqvar.VariantRecordWithCoordinate]=None,
            in_nodes:Set[PVGNode]=None, out_nodes:Set[PVGNode]=None,
            selenocysteins:List[seqvar.VariantRecordWithCoordinate]=None,
            cleavage:bool=False, truncated:bool=False, orf:List[int]=None,
            was_bridge:bool=False, pre_cleaved:bool=False, level:int=0,
            npop_collapsed:bool=False, cpop_collapsed:bool=False,
            upstream_indel_map:Dict[PVGNode, List[seqvar.VariantRecord]]=None,
            collapsed_variants:Dict[PVGNode, List[seqvar.VariantRecordWithCoordinate]]=None
            ):
        """ Construct a PVGNode object. """
        self.seq = seq
        self.variants = variants or []
        self.in_nodes = in_nodes or set()
        self.out_nodes = out_nodes or set()
        self.selenocysteins = selenocysteins or []
        self.cleavage = cleavage
        self.truncated = truncated
        self.orf = orf or [None, None]
        self.reading_frame_index = reading_frame_index
        self.was_bridge = was_bridge
        self.pre_cleave = pre_cleaved
        self.subgraph_id = subgraph_id
        self.level = level
        self.npop_collapsed = npop_collapsed
        self.cpop_collapsed = cpop_collapsed
        self.upstream_indel_map = upstream_indel_map or {}
        self.collapsed_variants = collapsed_variants or {}

    def __getitem__(self, index) -> PVGNode:
        """ get item """
        start, stop, _ = index.indices(len(self.seq))
        location = FeatureLocation(start=start, end=stop)
        seq = self.seq.__getitem__(index)
        variants = []

        for variant in self.variants:
            if variant.location.overlaps(location):
                variants.append(variant.shift(-start))

        npop_collapsed = self.npop_collapsed and start == 0
        cpop_collapsed = self.cpop_collapsed and stop == len(seq) or stop == -1

        if start == 0:
            upstream_indel_map = {k:copy.copy(v) for k,v in self.upstream_indel_map.items()}
        else:
            upstream_indel_map = []

        collapsed_variants = {}
        for upstream, variants in self.collapsed_variants.items():
            filtered_variants = []
            for v in variants:
                if v.location.overlaps(location):
                    v = v[start:stop]
                    v = v.shift(-start)
                    filtered_variants.append(v)
            if filtered_variants:
                collapsed_variants[upstream] = filtered_variants

        secs = []
        for sec in self.selenocysteins:
            if start <= sec.location.start < stop:
                secs.append(sec.shift(-start))

        return PVGNode(
            seq=seq,
            variants=variants,
            cleavage=self.cleavage,
            orf=self.orf,
            selenocysteins=secs,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            npop_collapsed=npop_collapsed,
            cpop_collapsed=cpop_collapsed,
            upstream_indel_map=upstream_indel_map,
            collapsed_variants=collapsed_variants
        )

    def add_out_edge(self, node:PVGNode) -> None:
        """ Add a outbound edge from this node.

        Args:
            node (PVGNode): The outboud node of the edge to add.
        """
        self.out_nodes.add(node)
        node.in_nodes.add(self)

    def remove_out_edge(self, node:PVGNode) -> None:
        """ Remove a outbound edge. It tries to remove the node, but won't
        raise any error if it doesn't exist.

        Args:
            node (PVGNode): The outbound node of the edge to remove.
        """
        self.out_nodes.discard(node)
        node.in_nodes.discard(self)

    def get_out_nodes(self) -> List[PVGNode]:
        """ Get outgoing nodes as a list. """
        return list(self.out_nodes)

    def get_in_nodes(self) -> List[PVGNode]:
        """ Get incoming nodes as a list """
        return list(self.in_nodes)

    def has_multiple_segments(self) -> bool:
        """ Whether the node has multiple segments, which is when the node
        is merged from several individual nodes. """
        n = len(self.seq.locations)
        for v in self.variants:
            if not (v.variant.is_fusion() \
                    or v.variant.is_circ_rna() \
                    or (v.variant.is_alternative_splicing() and not v.variant.is_deletion())):
                n += 1
        return n > 1

    def _get_nth_rf_index(self, i:int) -> int:
        """ """
        if (i > 0 or i < -1) and not self.has_multiple_segments():
            raise ValueError('Node does not have multiple segments')

        locations = [loc.query for loc in self.seq.locations]
        for v in self.variants:
            if not (v.variant.is_fusion() \
                    or v.variant.is_circ_rna() \
                    or (v.variant.is_alternative_splicing() and not v.variant.is_deletion())):
                locations.append(v.location)

        locations.sort()

        return locations[i].reading_frame_index

    def get_first_rf_index(self) -> int:
        """ Get the first fragment's reading frame index """
        return self._get_nth_rf_index(0)

    def get_second_rf_index(self) -> int:
        """ Get the second fragment's reading frame index """
        return self._get_nth_rf_index(1)

    def get_last_rf_index(self) -> int:
        """ Get the last fragment's reading frame index """
        return self._get_nth_rf_index(-1)

    def remove_out_edges(self) -> None:
        """ remove output nodes """
        for node in copy.copy(self.out_nodes):
            self.remove_out_edge(node)

    def is_inbond_of(self, node:PVGNode) -> bool:
        """ Checks if self is the inbond node of a given node. """
        return node in self.out_nodes

    def is_orphan(self) -> bool:
        """ Checks if the node is orphan (no inbond or outbond node) """
        return not self.in_nodes and not self.out_nodes

    def is_bridge(self) -> None:
        """ Check if this is a bridge node to another reading frame """
        for node in self.out_nodes:
            if not (not node.out_nodes and node.seq.seq == '*') and \
                    node.reading_frame_index != self.reading_frame_index:
                return True
        if self.was_bridge:
            return True
        return False

    def is_subgraph_end(self) -> bool:
        """ Check if this is the end node of a subgraph """
        return all(x.subgraph_id == self.subgraph_id for x in self.in_nodes) and \
            any(x.subgraph_id != self.subgraph_id for x in self.out_nodes)

    def is_subgraph_start(self) -> bool:
        """ Check if this is the start node of a subgraph """
        return any(x.subgraph_id != self.subgraph_id for x in self.in_nodes) and \
            all(x.subgraph_id == self.subgraph_id for x in self.out_nodes)

    def is_at_least_one_loop_downstream(self, other, subgraphs:SubgraphTree,
            circ_rna:circ.CircRNAModel) -> bool:
        """ Checks if is at least one loop downstream to the other node. """
        if self.seq.locations:
            i = self.seq.locations[-1].ref.start
            x_level = subgraphs[self.seq.locations[-1].ref.seqname].level
        else:
            for v in reversed(self.variants):
                if not v.variant.is_circ_rna():
                    i = math.floor(v.variant.location.start / 3)
                    x_level = subgraphs[v.location.seqname].level
                    break
            else:
                raise ValueError("Failed to find a non circRNA variant.")

        if other.seq.locations:
            j = other.seq.locations[0].ref.start
            y_level = subgraphs[other.seq.locations[0].ref.seqname].level
        else:
            for v in other.variants:
                if not v.variant.is_circ_rna():
                    j = math.floor(v.variant.location.start / 3)
                    y_level = subgraphs[v.location.seqname].level
                    break
            else:
                raise ValueError("Failed to find a non circRNA variant.")

        if x_level - y_level > 1:
            return True

        if x_level == y_level:
            return False

        if x_level - y_level == 1:
            for fragment in circ_rna.fragments:
                frag = FeatureLocation(
                    start=math.floor((fragment.location.start - 3) / 3),
                    end=math.floor(fragment.location.end / 3)
                )
                if i in frag and j in frag:
                    return i >= j
                if i in frag:
                    return False
                if j in frag:
                    return True
        raise ValueError('Locations not found from the fragments of the circRNA.')

    def has_any_in_bridge(self) -> None:
        """ Check if it has any incoming node that is bridge """
        return any(node.is_bridge() for node in self.in_nodes)

    def has_any_indel(self) -> None:
        """ Checks if there is any indel """
        return any(x.variant.is_indel() for x in self.variants)

    def get_variants_at(self, start:int, end:int=-1,
            upstream_cleavage_altering:bool=True,
            downstream_cleavage_altering:bool=True
            ) -> List[seqvar.VariantRecord]:
        """ Get the variant at position i """
        if end == -1:
            end = len(self.seq)
        if not -1 < start <= end <= len(self.seq):
            raise ValueError('start and end out of the range')
        variants = []
        location = FeatureLocation(start=start, end=end)
        for variant in self.variants:
            if not upstream_cleavage_altering and variant.upstream_cleavage_altering:
                continue
            if not downstream_cleavage_altering and variant.downstream_cleavage_altering:
                continue
            if variant.location.overlaps(location):
                variants.append(variant.variant)
        return variants

    def get_cleavage_gain_variants(self) -> List[seqvar.VariantRecord]:
        """ Get cleavage gain variants """
        cleavage_gain = []
        for variant in self.variants:
            if variant.location.end == len(self.seq) and not variant.is_silent:
                cleavage_gain.append(variant.variant)
        return cleavage_gain

    def get_cleavage_gain_from_downstream(self) -> List[seqvar.VariantRecord]:
        """ Get the variants that gains the cleavage by downstream nodes """
        cleavage_gain = []
        upstream_cleave_alts = [v.variant for v in self.variants
            if v.location.end == len(self.seq.seq)]
        for node in self.out_nodes:
            if not node.variants:
                return []
            if node.variants[0].location.start != 0:
                return []
            for v in node.variants:
                if v.location.start != 0:
                    break
                if not cleavage_gain and not v.is_silent \
                        and not v.downstream_cleavage_altering:
                    if v.upstream_cleavage_altering:
                        if not any(x == v.variant for x in upstream_cleave_alts):
                            continue
                    cleavage_gain.append(v.variant)
        return cleavage_gain

    def get_stop_lost_variants(self, stop_index:int) -> List[seqvar.VariantRecord]:
        """ Get stop lost variants """
        stop_lost = []
        stop_codon = FeatureLocation(start=stop_index, end=stop_index + 3)
        for variant in self.variants:
            if variant.variant.location.overlaps(stop_codon):
                stop_lost.append(variant.variant)
        return stop_lost

    def frames_shifted(self) -> int:
        """ Get total frames shifted of all variants """
        return sum([v.variant.frames_shifted() for v in self.variants]) % 3

    def find_reference_next(self) -> PVGNode:
        """ Find and return the next reference node. The next reference node
        is defined as the out node that has not variant, or not any variant
        that is not in this current node. """
        if not self.out_nodes:
            return None
        this_variants = [x.variant for x in self.variants]
        for node in self.out_nodes:
            if not node.variants:
                return node
            that_variants = [x.variant for x in node.variants]
            if not any(it not in this_variants for it in that_variants):
                return node
        raise ValueError('None of the out nodes is reference.')

    def find_reference_prev(self) -> PVGNode:
        """ Find and return the previous reference node. The previous reference
        node is defined as the inbond node that has not variant, or any variant
        that is not in this current node. """
        if not self.in_nodes:
            return None
        if len(self.in_nodes) == 1:
            return list(self.in_nodes)[0]
        this_variants = [x.variant for x in self.variants]
        for node in self.in_nodes:
            if not node.variants:
                return node
            that_variants = [x.variant for x in node.variants]
            if not any(it not in this_variants for it in that_variants):
                return node
        raise ValueError('None of the in nodes is reference.')

    def find_least_varaint_prev(self) -> PVGNode:
        """ Find and return the inbonding node that has the least number of
        variants """
        if not self.in_nodes:
            return None
        if len(self.in_nodes) == 1:
            return list(self.in_nodes)[0]

        def sort_func(this:PVGNode, that:PVGNode):
            x = this.variants
            y = that.variants
            if len(x) < len(y):
                return -1
            if len(x) > len(y):
                return 1
            if len(x) == 0 and len(y) == 0:
                return -1
            for i,j in zip(x,y):
                if i.location < j.location:
                    return -1
                if i.location > j.location:
                    return 1
            return -1

        return sorted(self.in_nodes, key=cmp_to_key(sort_func))[0]

    def find_least_variant_next(self) -> PVGNode:
        """ Find the outbond node that has the least number of variants """
        if not self.out_nodes:
            return None
        if len(self.out_nodes) == 1:
            return list(self.out_nodes)[0]

        try:
            return self.find_reference_next()
        except ValueError:
            pass

        self_variants = {v.variant for v in self.variants}
        def sort_func(this:PVGNode, that:PVGNode):
            x = [v for v in this.variants if v.variant not in self_variants]
            y = [v for v in that.variants if v.variant not in self_variants]
            if len(x) < len(y):
                return -1
            if len(x) > len(y):
                return 1
            if len(x) == 0 and len(y) == 0:
                return -1
            for i,j in zip(x,y):
                if i.location < j.location:
                    return -1
                if i.location > j.location:
                    return 1
            return -1

        return sorted(self.out_nodes, key=cmp_to_key(sort_func))[0]

    def split_selenocysteins(self, index:int):
        """ Split selenocysterines """
        left_secs:List[seqvar.VariantRecordWithCoordinate] = []
        right_secs:List[seqvar.VariantRecordWithCoordinate] = []
        for sec in self.selenocysteins:
            if sec.location.start < index:
                left_secs.append(sec)
            else:
                right_secs.append(sec.shift(-index))
        return left_secs, right_secs

    def split_node(self, index:int, cleavage:bool=False, pre_cleave:bool=False,
            pop_collapse:bool=False) -> PVGNode:
        """ Split the sequence at the given position, and create a new node
        as the outbound edge. Variants will also be adjusted. For example:

        KAPGCFP -> split at position 3 -> KAP-GCFP

        The node with GCFP is returned.

        Args:
            index (int): the position to split
            cleavage (bool): If true, the cleavage attribute of the right node
                is update to be True.

        Returns:
            The new created node with the right part of the original sequence.
        """
        left_seq = self.seq[:index]
        right_seq = self.seq[index:]

        left_variants = []
        right_variants = []
        for variant in self.variants:
            if variant.location.start < index:
                if variant.location.end <= index:
                    left_variants.append(variant)
                else:
                    left_variants.append(variant[:index])
            elif variant.location.start == index and cleavage \
                    and not variant.downstream_cleavage_altering:
                cleave_alts = variant[index:index+1]
                cleave_alts = cleave_alts.shift(-1)
                cleave_alts.downstream_cleavage_altering = True
                left_variants.append(cleave_alts)

            if variant.location.end > index:
                right_variants.append(variant.shift(-index))
            elif variant.location.end == index and cleavage \
                    and not variant.upstream_cleavage_altering:
                cleave_alts = variant.shift(-index)
                cleave_alts.upstream_cleavage_altering = True
                right_variants.append(cleave_alts)

        self.seq = left_seq
        self.variants = left_variants

        new_node = PVGNode(
            seq=right_seq,
            reading_frame_index=self.reading_frame_index,
            variants=right_variants,
            orf=self.orf,
            was_bridge=self.was_bridge,
            pre_cleaved=pre_cleave,
            subgraph_id=self.subgraph_id,
            level=self.level,
            cpop_collapsed=self.cpop_collapsed,
            truncated=self.truncated
        )
        new_node.orf = self.orf

        left_secs, right_secs = self.split_selenocysteins(index)
        self.selenocysteins = left_secs
        new_node.selenocysteins = right_secs

        if cleavage:
            new_node.cleavage = True

        if pop_collapse:
            self.cpop_collapsed = True
            new_node.npop_collapsed = True
            if new_node.is_bridge():
                self.was_bridge = True
        else:
            self.cpop_collapsed = False

        self.truncated = False

        while self.out_nodes:
            node = self.out_nodes.pop()
            self.remove_out_edge(node)
            new_node.add_out_edge(node)
        self.add_out_edge(new_node)
        return new_node

    def truncate_right(self, i:int) -> PVGNode:
        """ Truncate the right i nucleotides off. """
        right_seq = self.seq[i:]
        left_variants = []
        right_variants = []

        for variant in self.variants:
            if variant.location.start < i:
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))

        right_node = PVGNode(
            seq=right_seq,
            variants=right_variants,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            cpop_collapsed=self.cpop_collapsed
        )

        self.seq = self.seq[:i]
        self.variants = left_variants

        left_secs, right_secs = self.split_selenocysteins(i)
        self.selenocysteins = left_secs
        right_node.selenocysteins = right_secs

        return right_node

    def truncate_left(self, i:int) -> PVGNode:
        """ Truncate the left i nucleotides off. A new node with the left part
        of the sequences and variants associated is returned. The self node
        is updated with only the right part of the sequence and variants. """
        left_seq = self.seq[:i]
        left_variants = []
        right_variants = []

        for variant in self.variants:
            if variant.location.start < i:
                left_variants.append(variant)
            if variant.location.end > i:
                right_variants.append(variant.shift(-i))

        left_node = PVGNode(
            seq=left_seq,
            variants=left_variants,
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            npop_collapsed=self.npop_collapsed,
            upstream_indel_map=self.upstream_indel_map
        )
        self.upstream_indel_map = {}

        self.seq = self.seq[i:]
        self.variants = right_variants

        left_secs, right_secs = self.split_selenocysteins(i)
        left_node.selenocysteins = left_secs
        self.selenocysteins = right_secs

        return left_node

    def append_left(self, other:PVGNode) -> None:
        """ Combine the other node to the left. """
        self.seq = other.seq + self.seq
        variants = copy.copy(other.variants)
        for variant in self.variants:
            variants.append(variant.shift(len(other.seq.seq)))
        self.variants = variants
        self.upstream_indel_map = {k:copy.copy(v) for k,v in
            other.upstream_indel_map.items()}

        secs = copy.copy(other.selenocysteins)
        for sec in self.selenocysteins:
            sec = sec.shift(len(other.seq.seq))
            secs.append(sec)
        self.selenocysteins = secs

    def append_right(self, other:PVGNode) -> None:
        """ Combine the other node the the right. """
        new_seq = self.seq + other.seq
        variants = []
        for variant in self.variants:
            if not variant.downstream_cleavage_altering:
                variants.append(variant)
        self.variants = variants

        for variant in other.variants:
            self.variants.append(variant.shift(len(self.seq.seq)))

        self.seq = new_seq
        self.cpop_collapsed = other.cpop_collapsed
        self.truncated = other.truncated

        for sec in other.selenocysteins:
            sec = sec.shift(len(self.seq.seq))
            self.selenocysteins.append(sec)

    def find_start_index(self) -> int:
        """ Find the start amino acid position """
        return self.seq.seq.find('M')

    def copy(self, in_nodes:bool=True, out_nodes:bool=True) -> PVGNode:
        """ Create a copy of the node """
        new_in_nodes = copy.copy(self.in_nodes) if in_nodes else set()
        new_out_nodes = copy.copy(self.out_nodes) if out_nodes else set()
        return PVGNode(
            seq=self.seq,
            variants=copy.copy(self.variants),
            in_nodes=new_in_nodes,
            out_nodes=new_out_nodes,
            cleavage=self.cleavage,
            truncated=self.truncated,
            orf=self.orf,
            selenocysteins=copy.copy(self.selenocysteins),
            reading_frame_index=self.reading_frame_index,
            was_bridge=self.was_bridge,
            subgraph_id=self.subgraph_id,
            level=self.level,
            npop_collapsed=self.npop_collapsed,
            cpop_collapsed=self.cpop_collapsed,
            upstream_indel_map={k:copy.copy(v) for k,v in self.upstream_indel_map.items()}
        )

    def get_nearest_next_ref_index(self) -> int:
        """ get the nearest reference index """
        cur = self
        while not cur.seq.locations and cur.out_nodes:
            cur = cur.find_reference_next()

        if not cur.seq.locations:
            return -1
        return int(cur.seq.locations[0].ref.start)

    def get_orf_start(self, i:int=0) -> int:
        """ get the ORF start position at the reference transcript sequence
        given M is found at the query postion i.
        """
        k = self.reading_frame_index
        if self.seq.locations:
            for loc in self.seq.locations:
                if i < loc.query.start:
                    return loc.ref.start * 3 + k
                if loc.query.start <= i < loc.query.end:
                    return (i - loc.query.start + loc.ref.start) * 3 + k

        out_node = self.find_least_variant_next()
        if str(out_node.seq.seq) == '*' and not out_node.out_nodes:
            return -1
        if out_node.seq.locations:
            return out_node.seq.locations[0].ref.start * 3 + k
        # raise ValueError('Can not find ORF')
        return -1

    def get_subgraph_id_at(self, i:int) -> str:
        """ Get the subgraph ID at the given position """
        for loc in self.seq.locations:
            if i in loc.query:
                return loc.ref.seqname
            if loc.query.start > i:
                break
        for v in self.variants:
            if i in v.location:
                return v.location.seqname
            if v.location.start > i:
                break
        raise ValueError('Failed to find the subgraph ID.')

    def is_identical(self, other:PVGNode) -> bool:
        """ Checks if two nodes have the same sequence and same outbond node """
        return self.seq == other.seq and \
            self.out_nodes == other.out_nodes and \
            self.cleavage == other.cleavage and \
            self.reading_frame_index == other.reading_frame_index

    def is_less_mutated_than(self, other:PVGNode) -> bool:
        """ Checks if this node has less mutation than the other """
        if self.level != other.level:
            return self.level < other.level

        if len(self.variants) != len(other.variants):
            return len(self.variants) < len(other.variants)

        for x, y in zip(self.variants, other.variants):
            if x.location != y.location:
                return x.location < y.location

        for x, y in zip(self.variants, other.variants):
            # SNV is higher than RES
            if x.variant.type == 'SNV' and y.variant.type == 'RNAEditingSite':
                return True
            if x.variant.type == 'RNAEditingSite' and y.variant.type == 'SNV':
                return False
            # Otherwise don't really care about the order, just to make sure
            # the result is reproducible
            if x.variant.type != y.variant.type:
                return x.variant.type > y.variant.type

            if x.variant.alt != y.variant.type:
                return x.variant.alt > y.variant.type
        return True

    def transfer_in_nodes_to(self, other:PVGNode):
        """ transfer all in nodes of current node to the other """
        for in_node in copy.copy(self.in_nodes):
            if in_node in other.in_nodes:
                continue
            in_node.add_out_edge(other)

    def is_already_cleaved(self):
        """ Check if the node is already cleaved """
        return (self.cleavage and not self.pre_cleave) and \
            all(x.cleavage and not x.pre_cleave for x in self.out_nodes) and\
            all(all(y.cleavage and not y.pre_cleave for y in x.in_nodes)
                for x in self.out_nodes)

    def get_downstream_stop_altering_variants(self) -> Set[VariantRecord]:
        """ Get downstream stop altering variants """
        final_variants = set()
        for out_node in self.out_nodes:
            if out_node.seq.seq == '*':
                stop_alts = set()
                stop_alts.update([x.variant for x in out_node.variants
                    if x.is_stop_altering and not x.downstream_cleavage_altering])
                stop_alts.update(out_node.upstream_indel_map.get(self, []))
                if not stop_alts:
                    return set()
                if not final_variants:
                    final_variants = stop_alts
        return final_variants

    def has_variant_at(self, start:int, end:int) -> bool:
        """ Checks if the node has any variant at a given position """
        loc = FeatureLocation(start=start, end=end)
        for variant in self.variants:
            if variant.location.overlaps(loc):
                return True
        return False

    def is_missing_variant(self, variant:seqvar.VariantRecord, upstream:PVGNode=None
            ) -> bool:
        """ Checks if in the coordinates of the node if any variant from a
        given set is missing. This is used for circRNA to ensure that the
        same variant is either on or off in each loop. """
        locs = []
        if self.seq.locations:
            for loc in self.seq.locations:
                tx_loc = FeatureLocation(
                    start=loc.get_ref_dna_start(),
                    end=loc.get_ref_dna_end()
                )
                locs.append(tx_loc)
        if self.variants:
            for v in self.variants:
                if v.variant.is_circ_rna():
                    continue
                locs.append(v.variant.location)
        variants = {x.variant for x in self.variants}
        if upstream:
            upstream_indels = self.upstream_indel_map.get(upstream, [])
            variants.update(v for v in upstream_indels)
        for loc in locs:
            if loc.overlaps(variant.location) \
                    and not any(variant == v for v in variants):
                return True
        return False

    def is_missing_any_variant(self, variants:Iterable[seqvar.VariantRecord],
            upstream:PVGNode=None) -> bool:
        """ Checks if the node is missing any of the variants. """
        return any(self.is_missing_variant(v, upstream) for v in variants)

    def has_variants_not_in(self, variants:Set[seqvar.VariantRecord]) -> bool:
        """ Checks if the node has any variants that is not in a given set. """
        if not variants:
            return False
        var_start = sorted(variants, key=lambda x: x.location.start)[0].location.start
        var_end = sorted(variants, key=lambda x: x.location.end)[-1].location.end
        var_range = FeatureLocation(start=var_start, end=var_end)
        for v in self.variants:
            if v.variant.is_circ_rna():
                continue
            if v.variant.location.overlaps(var_range) \
                    and v.variant not in variants:
                return True
        return False

    def is_hybrid_node(self, tree:SubgraphTree) -> bool:
        """ Checks if the node is hybrid. A hybrid node is when two parts of
        the node sequence are from different subgraphs, and have the different
        states. """
        variants = {v.variant for v in self.variants}

        seq_iter = iter(self.seq.locations)
        var_iter = iter(self.variants)

        cur_seq = next(seq_iter, None)
        cur_var = next(var_iter, None)

        start, end = 0, 0
        cur_level = -1

        while cur_seq or cur_var:
            if cur_var and cur_var.variant.is_circ_rna():
                cur_var = next(var_iter, None)
                continue

            seq_level, var_level = -1, -1
            seq_start, seq_end = -1, -1
            var_start, var_end = -1, -1

            if cur_seq:
                seq_level = tree[cur_seq.ref.seqname].level
                seq_start = cur_seq.query.start
                seq_end = cur_seq.query.end

            if cur_var:
                var_level = tree[cur_var.location.seqname].level
                var_start = cur_var.location.start
                var_end = cur_var.location.end

            seq_is_smaller = (not cur_var ) \
                or (cur_seq and (seq_level < var_level or seq_start < var_start))

            if seq_is_smaller:
                if cur_level == -1:
                    cur_level = seq_level
                    start, end = seq_start, seq_end
                elif seq_level > cur_level:
                    seg = self[start:end]
                    if seg.is_missing_any_variant(variants):
                        return True
                    cur_level = seq_level
                    start, end = seq_start, seq_end
                else:
                    end = seq_end

                if cur_seq is self.seq.locations[-1] and not cur_var:
                    seg = self[start:seq_end]
                    return seg.is_missing_any_variant(variants)

                cur_seq = next(seq_iter, None)

            else:
                if cur_level == -1:
                    cur_level = var_level
                    start, end = var_start, var_end
                elif var_level > cur_level:
                    seg = self[start:end]
                    if seg.is_missing_any_variant(variants):
                        return True
                    cur_level = var_level
                    start, end = var_start, var_end
                else:
                    end = var_end

                if cur_var is self.variants[-1] and not cur_seq:
                    seg = self[start:var_end]
                    return seg.is_missing_any_variant(variants)

                cur_var = next(var_iter, None)

        return False

    def any_unaccounted_downstream_cleavage_or_stop_altering(self,
            variants:Iterable[seqvar.VariantRecord]) -> bool:
        """ Checks if the node has any unaccounted cleavage or stop altering
        mutation from a given list. The first amino acid from the downstream
        node and the last amino acid from the upstream one will be checked if
        they contain any variant in the list of variant provided and not having
        any variant that are not in the list.
        """
        variants = {v for v in variants if not v.is_circ_rna()}

        boundary_nodes:List[PVGNode] = []

        # stop and downstream cleavage gain
        if len(self.get_out_nodes()) == 1:
            downstream = self.get_out_nodes()[0]
            if not (downstream.seq.seq == '*' and not downstream.get_out_nodes()):
                boundary_nodes.append(downstream[:1])

        # upstream cleavage gain cleavage gain
        if len(self.get_in_nodes()) == 1:
            upstream = self.get_in_nodes()[0]
            if upstream.seq is not None:
                boundary_nodes.append(upstream[-1:])

        for node in boundary_nodes:
            if node.is_missing_any_variant(variants, self):
                return True
            if node.has_variants_not_in(variants):
                return True
        return False

    def fix_selenocysteines(self,
            sect_variants:List[seqvar.VariantRecordWithCoordinate]) -> None:
        """ Fix Selenocysteines from the amin acid sequence. """
        iter_loc = iter(self.seq.locations)
        iter_sec = iter(sect_variants)
        loc = next(iter_loc, None)
        sect = next(iter_sec, None)

        sects:List[seqvar.VariantRecordWithCoordinate] = []

        while loc and sect:
            rf_index = loc.query.reading_frame_index
            if rf_index != sect.location.start % 3:
                sect = next(iter_sec, None)
                continue

            dna_start = loc.get_ref_dna_start()
            dna_end = loc.get_ref_dna_end()

            dna_loc = FeatureLocation(dna_start, dna_end)
            if dna_loc.is_superset(sect.location):
                k = loc.query.start + int((sect.location.start - dna_loc.start) / 3)
                sect_local = seqvar.VariantRecordWithCoordinate(
                    location=FeatureLocation(k, k+1),
                    variant=sect.variant
                )
                if self.seq.seq[k] != '*':
                    err.warning(
                        'The codon at the given Selenocysteine position is not' +
                        ' a stop codon.'
                    )
                sects.append(sect_local)

                sect = next(iter_sec, None)
                loc = next(iter_loc, None)
                continue

            if dna_loc > sect.location:
                sect = next(iter_sec, None)
            else:
                loc = next(iter_loc, None)

        if not sects:
            return

        for i,sect in enumerate(sects):
            k = sect.location.start
            if i == 0:
                new_seq = self.seq.seq[:k]
            elif sects[i-1] + 1 < k:
                new_seq += self.seq.seq[sects[i-1] + 1:k]
            new_seq += 'U'
            if k == sects[-1]:
                if k + 1 < len(self.seq.seq):
                    new_seq += self.seq.seq[k+1:]

        self.selenocysteins = sects
        self.seq.seq = new_seq
