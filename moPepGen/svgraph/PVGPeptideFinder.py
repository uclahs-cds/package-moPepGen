""" Module for variant peptide dict """
from __future__ import annotations
from collections import deque
import itertools
import copy
from typing import Deque, Dict, Iterable, List, Set, Tuple, TYPE_CHECKING
from Bio.Seq import Seq
from Bio import SeqUtils
from moPepGen import aa, get_equivalent, seqvar
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.svgraph.PVGOrf import PVGOrf
from moPepGen.aa import VariantPeptideIdentifier as vpi
from moPepGen.svgraph.SubgraphTree import SubgraphTree
from moPepGen.seqvar import create_variant_w2f


if TYPE_CHECKING:
    from moPepGen.seqvar.VariantRecord import VariantRecord
    from moPepGen.params import CleavageParams
    from moPepGen.circ import CircRNAModel

class PVGPeptideMetadata():
    """ Variant peptide metadata """
    def __init__(self, label:str=None, orf:Tuple[int,int]=None,
            is_pure_circ_rna:bool=False, has_variants:bool=False,
            segments:List[PeptideSegment]=None, check_orf:bool=False):
        """  """
        self.label = label
        self.orf = orf
        self.is_pure_circ_rna = is_pure_circ_rna
        self.has_variants = has_variants
        self.segments = segments or []
        self.check_orf = check_orf

    def get_key(self) -> str:
        """ get key """
        return f"{self.label}|{self.orf}" if self.check_orf else self.label

class PeptideSegment:
    """ Peptide segment """
    def __init__(self, query:FeatureLocation, ref:FeatureLocation, feature_type:str,
            feature_id:str, variant_id:str):
        """ constructor """
        self.query = query
        self.ref = ref
        self.feature_type = feature_type
        self.feature_id = feature_id
        self.variant_id = variant_id

    def __repr__(self):
        """ str """
        return f"<PeptideSegment query={self.query} ref={self.ref}" +\
            f" feature_type={self.feature_type} feature_id={self.feature_id}" +\
            f" variant_id={self.variant_id}>"

    def merge(self, other:PeptideSegment) -> PeptideSegment:
        """ merge """
        query = FeatureLocation(
            start=self.query.start,
            end=other.query.end,
            start_offset=self.query.start_offset,
            end_offset=other.query.end_offset
        )
        if self.ref:
            ref = FeatureLocation(
                start=self.ref.start,
                end=other.ref.end,
                start_offset=self.ref.start_offset,
                end_offset=other.ref.end_offset
            )
        else:
            ref = None
        return PeptideSegment(
            query=query,
            ref=ref,
            feature_type=self.feature_type,
            feature_id=self.feature_id,
            variant_id=self.variant_id
        )

    def is_adjacent(self, other:PeptideSegment) -> bool:
        """ is adjacent """
        if (self.variant_id is None) != (other.variant_id is None) \
                or (self.variant_id is not None \
                    and self.variant_id != other.variant_id) \
                or self.query.end != other.query.start \
                or self.feature_type != other.feature_type \
                or self.feature_id != other.feature_id:
            return False
        if self.ref:
            this_end = self.ref.end * 3 + self.ref.end_offset - self.query.end_offset
            that_start = other.ref.start * 3 + other.ref.start_offset + other.query.start_offset
            if this_end != that_start:
                return False
        return True

    def to_line(self):
        """ to line """
        fields = [
            str(self.query.start),
            str(self.query.end),
            self.feature_type or '.',
            self.feature_id or '.'
        ]
        if self.ref:
            ref_start = self.ref.start * 3 + self.ref.start_offset
            ref_end = self.ref.end * 3 + self.ref.end_offset
            fields += [
                str(ref_start),
                str(ref_end)
            ]
        else:
            fields += ['.', '.']
        fields += [
            str(self.query.start_offset),
            str(self.query.end_offset)
        ]
        if self.variant_id:
            fields.append(self.variant_id)
        else:
            fields.append('.')
        return '\t'.join(fields)

class AnnotatedPeptideLabel:
    """ Annotated peptide label """
    def __init__(self, label:str, segments:List[PeptideSegment]):
        """ constructor """
        self.label = label
        self.segments = segments

    def to_lines(self) -> str:
        """ to line """
        return [f"{self.label}\t{seg.to_line()}" for seg in self.segments]

class PVGNodePath():
    """ Helper class when calling for miscleavage peptides. The nodes contained
    by this class can be joined to make a miscleavage peptide. """
    def __init__(self, nodes:List[PVGNode], additional_variants:Set[VariantRecord],
            global_variants:Set[VariantRecord]=None):
        """ Constructor """
        self.nodes = nodes
        self.additional_variants = additional_variants
        self.global_variants = global_variants or set()
        self._len = None
        self._variants:Dict[str, VariantRecord] = None

    def __len__(self) -> int:
        """ Get length of node sequences """
        if self._len is None:
            self._len = sum(len(x.seq.seq) for x in self.nodes)
        return self._len

    def __copy__(self) -> PVGNodePath:
        """ copy """
        new_path = PVGNodePath(
            nodes=copy.copy(self.nodes),
            additional_variants=copy.copy(self.additional_variants)
        )
        new_path._len = self._len
        return new_path

    def is_too_short(self, param:CleavageParams) -> bool:
        """ Checks whether the sequence is too short """
        return len(self) < param.min_length

    def is_too_long(self, param:CleavageParams) -> bool:
        """ Checks whether the sequence is too long """
        return not(
            len(self) <= param.max_length
            or (
                self.nodes[0].seq.seq.startswith('M')
                and len(self) <= param.max_length + 1
            )
        )

    def has_trailing_selenocysteins(self):
        """ Checks whether the seires of nodes has any selenocystein in any of
        the trailing non C-pop-collapsed node. """
        for node in reversed(self.nodes):
            if node is not self.nodes[-1] and not node.cpop_collapsed:
                return False
            if node.selenocysteines:
                return True
        return False

    @property
    def variants(self):
        """ Gather variants from node path """
        if self._variants is None:
            self._variants = self.gather_variants()
        return self._variants

    def gather_variants(self) -> Dict[str, VariantRecord]:
        """ Gather variants """
        variants:Dict[str, VariantRecord] = {}
        for n in self.nodes:
            for v in n.variants:
                if v.variant.id not in variants:
                    variants[v.variant.id] = v.variant
        for v in self.additional_variants:
            if v.id not in variants:
                variants[v.id] = v
        return variants

    def append(self, node:PVGNode) -> None:
        """ Append node """
        assert self.nodes[-1].is_inbond_of(node)
        self.nodes.append(node)
        if self._len is not None:
            self._len += len(node.seq.seq)

    def add_additional_variants(self, variants:Iterable[VariantRecord]):
        """ Add additional variants """
        self.additional_variants.update(variants)
        if self._variants is not None:
            for v in variants:
                if v.id not in self._variants:
                    self._variants[v.id] = v

    def is_subgraph_spanning(self) -> bool:
        """ Checks if the nodes are spanning over subgraphs. """
        subgraph_ids = set().union(*[n.get_subgraph_id_set() for n in self.nodes])
        return len(subgraph_ids) > 1

    def number_of_nodes(self) -> int:
        """ len """
        return len(self.nodes)

    def has_any_potential_island_variant(self, node:PVGNode) -> bool:
        """ Checks if the given node has any non global variant """
        for v in node.variants:
            if v.variant != node.global_variant and v.variant not in self.global_variants:
                return True
        return False

class PVGNodePathArchipel(PVGNodePath):
    """ PVGNodePathArchipel """
    def __init__(self, nodes:List[PVGNode], additional_variants:Set[VariantRecord],
        nflanking:List[PVGNode]=None, cflanking:List[PVGNode]=None, flanking_size:int=9,
        islands:List[int]=None, global_variants:Set[VariantRecord]=None):
        super().__init__(
            nodes=nodes,
            additional_variants=additional_variants,
            global_variants=global_variants
        )
        self.nflanking = nflanking or []
        self.cflanking = cflanking or []
        self.flanking_size = flanking_size
        self.islands = islands or []

    def is_nflanking_full(self):
        """ Is nflanking full """
        return len(self.nflanking) >= self.flanking_size

    def is_cflanking_full(self):
        """ Is cflanking full """
        return len(self.cflanking >= self.flanking_size)

    def has_any_island(self):
        """ Has any island """
        return len(self.islands) > 0

    def __copy__(self) -> PVGNodePathArchipel:
        """ copy """
        new_path:PVGNodePathArchipel = super().__copy__()
        new_path.__class__ = self.__class__
        new_path.nflanking = copy.copy(self.nflanking)
        new_path.cflanking = copy.copy(self.cflanking)
        new_path.flanking_size = self.flanking_size
        new_path.islands = copy.copy(self.islands)
        new_path.global_variants = copy.copy(self.global_variants)
        return new_path

    def to_reef(self) -> PVGNodePathReef:
        """ Convert to PVGNodePathReef """
        return PVGNodePathReef(
            nodes=self.nodes,
            additional_variants=self.additional_variants,
            global_variants=self.global_variants,
            flanking_size=self.flanking_size
        )

    def is_valid_path(self) -> bool:
        """ Check if path is valid """
        return self.has_any_island()

class PVGNodePathReefMetadata:
    """ PVGNodePathReefMetadata """
    def __init__(self, subgraph_id:int, global_variant:VariantRecord=None,
            nodes:List[int]=None):
        """ constructor """
        self.subgraph_id = subgraph_id
        self.global_variant = global_variant
        self.nodes = nodes or []

def append_reef(reefs:List[PVGNodePathReefMetadata], i:int, node:PVGNode):
    """ append reef """
    if not reefs:
        reef = PVGNodePathReefMetadata(
            subgraph_id=node.subgraph_id,
            global_variant=node.global_variant,
            nodes=[i]
        )
        reefs.append(reef)
    elif node.subgraph_id == reefs[-1].subgraph_id:
        reefs[-1].nodes.append(i)
    else:
        reef = PVGNodePathReefMetadata(
            subgraph_id=node.subgraph_id,
            global_variant=node.global_variant,
            nodes=[i]
        )
        reefs.append(reef)
    return reefs

class PVGNodePathReef(PVGNodePath):
    """ PVGNodePathReef """
    def __init__(self, nodes:List[PVGNode], additional_variants:Set[VariantRecord],
            reefs:List[PVGNodePathReefMetadata]=None, global_variants:Set[VariantRecord]=None,
            flanking_size:int=9):
        """ Constructor """
        super().__init__(
            nodes=nodes,
            additional_variants=additional_variants,
            global_variants=global_variants
        )
        # (subgraph_id, nodes indices)
        self.reefs = reefs or self.find_reefs()
        self.flanking_size = flanking_size

    def find_reefs(self) -> List[PVGNodePathReefMetadata]:
        """ find all reefs """
        reefs:List[PVGNodePathReefMetadata] = []
        for i, node in enumerate(self.nodes):
            append_reef(reefs, i, node)
        return reefs

    def is_last_reef_full(self) -> bool:
        """ Check if the last reef is full """
        return self.reefs[-1].global_variant is None \
            and len(self.reefs[-1].nodes) >= self.flanking_size

    def add_node(self, node:PVGNode):
        """ Add node """
        self.nodes.append(node)
        if node.global_variant:
            self.global_variants.add(node.global_variant)
        i = len(self.nodes)
        append_reef(self.reefs, i, node)

    def __copy__(self) -> PVGNodePathReef:
        """ copy """
        new_path:PVGNodePathReef = super().__copy__()
        new_path.__class__ = self.__class__
        new_path.reefs = copy.deepcopy(self.reefs)
        new_path.flanking_size = self.flanking_size
        new_path.global_variants = copy.copy(self.global_variants)
        return new_path

    def is_valid_path(self) -> bool:
        """ Check if path is valid """
        return True

TypeVariantPeptideMetadataMap = Dict[Seq, Dict[str, PVGPeptideMetadata]]

class PVGCandidateNodePaths():
    """ Helper class for looking for peptides with miscleavages. This class
    defines the collection of nodes in sequence that starts from the same node,
    with up to X allowed miscleavages that makes the miscleaved peptide
    sequences.

    Attributes:
        - `data` (Deque[List[PVGNode]]): The collection of nodes that make all
          miscleaved sequences.
        - `cleavage_params` (CleavageParams): Cleavage related parameters.
        - `orf` (Tuple[int,int]): The ORF start and end locations.
        - `tx_id` (str): Transcript ID.
        - `leading_node` (PVGNode): The start node that the miscleaved peptides
          are called from. This node must present in the PVG graph.
    """
    def __init__(self, data:Deque[PVGNodePath],
            cleavage_params:CleavageParams, orfs:List[PVGOrf]=None, tx_id:str=None,
            gene_id:str=None, leading_node:PVGNode=None, subgraphs:SubgraphTree=None,
            is_circ_rna:bool=False):
        """ constructor """
        self.data = data
        self.orfs = orfs or []
        self.tx_id = tx_id
        self.gene_id = gene_id
        self.cleavage_params = cleavage_params
        self.leading_node = leading_node
        self.subgraphs = subgraphs
        self.is_circ_rna = is_circ_rna

    def is_valid_seq(self, seq:Seq, pool:Set[Seq], denylist:Set[str]) -> bool:
        """ Check whether the seq is valid """
        if seq in pool:
            return True
        min_mw = self.cleavage_params.min_mw
        return self.seq_has_valid_size(seq) \
            and seq not in denylist \
            and 'X' not in seq \
            and SeqUtils.molecular_weight(seq, 'protein') >= min_mw

    def seq_has_valid_size(self, seq:Seq=None, size:int=None) -> bool:
        """ Check whether the seq has valid length """
        if seq is None and size is None:
            raise ValueError
        min_length = self.cleavage_params.min_length
        max_length = self.cleavage_params.max_length
        if size is None:
            size = len(seq)

        return min_length <= size <= max_length

    def create_peptide_segments(self, nodes:List[PVGNode]) -> List[PeptideSegment]:
        """ create peptide segments """
        ref_segments:List[PeptideSegment] = []
        var_segments:List[PeptideSegment] = []
        offset = 0
        branch_variants = set()
        for node in nodes:
            for loc in node.seq.locations:
                branch = self.subgraphs[loc.ref.seqname]
                if branch.variant:
                    branch_variants.add(branch.variant.id)
                seg = PeptideSegment(
                    query=loc.query.shift(offset),
                    ref=loc.ref,
                    feature_type=branch.feature_type,
                    feature_id=branch.feature_id,
                    variant_id=branch.variant.id if branch.variant else None
                )
                if ref_segments and ref_segments[-1].is_adjacent(seg):
                    ref_segments[-1] = ref_segments[-1].merge(seg)
                    continue
                ref_segments.append(seg)
            for variant in node.variants:
                if variant.variant.id in branch_variants:
                    continue
                if variant.upstream_cleavage_altering:
                    continue
                if variant.downstream_cleavage_altering:
                    continue
                if len(variant.location) < 1:
                    continue
                seg = PeptideSegment(
                    query=variant.location.shift(offset),
                    ref=None,
                    feature_type=None,
                    feature_id=None,
                    variant_id=variant.variant.id
                )
                if var_segments and var_segments[-1].is_adjacent(seg):
                    var_segments[-1] = var_segments[-1].merge(seg)
                    continue
                var_segments.append(seg)
            offset += len(node.seq.seq)
        segments = ref_segments + var_segments
        segments.sort(key=lambda x:x.query)
        return segments

    def join_peptides(self, pool:TypeVariantPeptideMetadataMap,
            check_variants:bool, additional_variants:List[VariantRecord],
            denylist:Set[str], is_start_codon:bool=False,
            circ_rna:CircRNAModel=None, truncate_sec:bool=False,
            check_external_variants:bool=True, check_orf:bool=False
            ) -> Iterable[Tuple[Seq, PVGPeptideMetadata]]:
        """ join miscleaved peptides and update the peptide pool.

        Args:
            - `check_variants` (bool): When true, only peptides that carries at
              least 1 variant are kept. And when false, all unique peptides are
              reported (e.g. noncoding).
            - `additional_variants` (List[VariantRecord]): Additional variants,
              e.g., start gain, cleavage gain, stop lost.
            - `denylist` (Set[str]): Peptide sequences that should be excluded.
            - `is_start_codon` (bool): Whether the node contains start codon.
        """
        for series in self.data:
            queue = series.nodes
            metadata = PVGPeptideMetadata(check_orf=check_orf)
            seqs_to_join:List[Seq] = []
            size:int = 0
            variants:Dict[str,VariantRecord] = {}
            in_seq_variants:Dict[str,VariantRecord] = {}

            selenocysteines = []

            for i, node in enumerate(queue):
                other = str(node.seq.seq)
                if truncate_sec:
                    if size > 0:
                        selenocysteines += [x.shift(size) for x in node.selenocysteines]
                    else:
                        selenocysteines +=  node.selenocysteines
                seqs_to_join.append(other)
                size += len(other)

                if check_variants:
                    for variant in node.variants:
                        if variant.is_silent:
                            continue
                        if i > 0 and variant.upstream_cleavage_altering:
                            continue
                        if variant.variant.id not in variants:
                            variants[variant.variant.id] = variant.variant
                        if not variant.variant.is_circ_rna():
                            if variant.variant.id not in in_seq_variants:
                                in_seq_variants[variant.variant.id] = variant.variant
                    if i < len(queue) - 1:
                        _node = self.leading_node if i == 0 else node
                        indels = queue[i + 1].upstream_indel_map.get(_node)
                        if indels:
                            for variant in indels:
                                if variant.id not in variants:
                                    variants[variant.id] = variant

            valid_orf = None
            for orf in self.orfs:
                if self.is_circ_rna \
                        and (any(v.is_circ_rna() for v in variants.values()) \
                        or any(v.is_circ_rna() for v in orf.start_gain)):
                    if orf.is_valid_orf_to_node_path(queue, self.subgraphs, circ_rna):
                        if any(n.is_missing_any_variant(in_seq_variants.values()) for n in queue):
                            continue
                        metadata.orf = tuple(orf.orf)
                        valid_orf = orf
                        break
                else:
                    metadata.orf = tuple(orf.orf)
                    valid_orf = orf
                    break

            if valid_orf is None:
                continue

            cleavage_gain_down = queue[-1].get_cleavage_gain_from_downstream()
            for v in cleavage_gain_down:
                if v not in variants:
                    variants[v.id] = v

            if check_variants:
                for v in valid_orf.start_gain:
                    if v.id not in variants:
                        variants[v.id] = v
                for v in additional_variants:
                    if v.id not in variants:
                        variants[v.id] = v
                for v in valid_orf.cleavage_gain:
                    if v.id not in variants:
                        variants[v.id] = v
                for v in series.additional_variants:
                    if v.id not in variants:
                        variants[v.id] = v

            if self.is_circ_rna \
                    and any(v.is_circ_rna() for v in variants.values()) \
                    and queue[-1].any_unaccounted_downstream_cleavage_or_stop_altering(
                        {x for x in variants.values() if not x.is_circ_rna()}):
                continue

            if size == 0:
                continue

            if check_variants:
                if check_external_variants:
                    if not variants:
                        continue
                else:
                    if not (variants or selenocysteines):
                        continue

            if not selenocysteines \
                    and not self.seq_has_valid_size(size=size) \
                    and not (seqs_to_join[0].startswith('M')
                                and self.seq_has_valid_size(size=size-1)):
                continue

            seq = Seq(''.join(seqs_to_join))
            is_in_denylist = seq in denylist and (not is_start_codon or seq[1:] in denylist)
            if not seq in pool and is_in_denylist:
                continue

            metadata.is_pure_circ_rna = self.is_circ_rna \
                and len(variants) == 1 \
                and list(variants.values())[0].is_circ_rna()

            seqs = self.translational_modification(seq, metadata, denylist,
                variants.values(), is_start_codon, selenocysteines,
                check_variants, check_external_variants, pool, queue
            )
            for seq, metadata in seqs:
                yield seq, metadata


    def translational_modification(self, seq:Seq, metadata:PVGPeptideMetadata,
            denylist:Set[str], variants:Set[VariantRecord], is_start_codon:bool,
            selenocysteines:List[seqvar.VariantRecordWithCoordinate],
            check_variants:bool, check_external_variants:bool, pool:Set[Seq],
            nodes:List[PVGNode]
            ) -> Iterable[Tuple[Seq,PVGPeptideMetadata]]:
        """ Apply any modification that could happen during translation. The
        kinds of modifications that could happen are:
        1. Leading Methionine truncation.
        2. Selenocysteine truncation.
        """
        if variants or not check_variants:
            is_valid = self.is_valid_seq(seq, pool, denylist)

            is_valid_start =  is_start_codon and seq.startswith('M') and\
                self.is_valid_seq(seq[1:], pool, denylist)

            if is_valid or is_valid_start:
                cur_metadata = copy.copy(metadata)
                label = vpi.create_variant_peptide_id(
                    transcript_id=self.tx_id, variants=variants, orf_id=None,
                    gene_id=self.gene_id
                )
                cur_metadata.label = label
                cur_metadata.has_variants = bool(variants)

                if is_valid:
                    cur_metadata_2 = copy.copy(cur_metadata)
                    cur_metadata_2.segments = self.create_peptide_segments(nodes)
                    yield seq, cur_metadata_2

                if is_valid_start:
                    cur_seq=seq[1:]
                    cur_nodes = list(nodes)
                    cur_nodes[0] = cur_nodes[0].copy()
                    cur_nodes[0].truncate_left(1)
                    cur_metadata.segments = self.create_peptide_segments(cur_nodes)
                    yield cur_seq, cur_metadata

        # Selenocysteine termination
        for sec in selenocysteines:
            seq_mod = seq[:sec.location.start]
            is_valid = self.is_valid_seq(seq_mod, pool, denylist)
            is_valid_start = is_start_codon and seq_mod.startswith('M') and\
                self.is_valid_seq(seq_mod[1:], pool, denylist)

            if is_valid or is_valid_start:
                cur_metadata = copy.copy(metadata)
                cur_variants = [v for v in variants if v.location.end
                    <= sec.variant.location.start]
                if check_variants and check_external_variants and not cur_variants:
                    continue
                cur_variants.append(sec.variant)
                label = vpi.create_variant_peptide_id(
                    transcript_id=self.tx_id, variants=cur_variants, orf_id=None,
                    gene_id=self.gene_id
                )
                cur_metadata.label = label
                cur_metadata.has_variants = bool(cur_variants)

                if is_valid:
                    cur_metadata_2 = copy.copy(cur_metadata)
                    cur_nodes = []
                    cut_offset = sec.location.start
                    for node in nodes:
                        if cut_offset == 0:
                            cur_nodes.append(node)
                            continue
                        if len(node.seq.seq) > cut_offset:
                            node = node.copy()
                            node.truncate_left(cut_offset)
                            cur_nodes.append(node)
                            cut_offset = 0
                        else:
                            cut_offset = max(0, cut_offset - len(node.seq.seq))
                            continue
                    cur_metadata_2.segments = self.create_peptide_segments(cur_nodes)
                    yield seq_mod, cur_metadata_2

                if is_valid_start:
                    cur_seq=seq_mod[1:]
                    cur_nodes[0] = cur_nodes[0].copy()
                    cur_nodes[0].truncate_left(1)
                    cur_metadata.segments = self.create_peptide_segments(cur_nodes)
                    yield cur_seq, cur_metadata

    @staticmethod
    def any_overlaps_and_all_missing_variant(nodes:Iterable[PVGNode], variant:VariantRecord):
        """ Whether any node overlaps with a given variant and all nodes are
        missing it. Should be only used for circRNA. """
        any_overlap = False
        all_missing = True
        for node in nodes:
            locs = []
            if node.seq.locations:
                for loc in node.seq.locations:
                    tx_loc = FeatureLocation(
                        start=int(loc.ref.start) * 3 + node.reading_frame_index,
                        end=int(loc.ref.end) * 3 + node.reading_frame_index
                    )
                    locs.append(tx_loc)
            if node.variants:
                locs += [x.variant.location for x in node.variants
                    if not x.variant.is_circ_rna()]
            node_variants = {x.variant for x in node.variants}
            if not any(loc.overlaps(variant.location) for loc in locs):
                continue
            any_overlap = True
            for loc in locs:
                if loc.overlaps(variant.location) \
                        and any(variant == v for v in node_variants):
                    all_missing = False
        return any_overlap and all_missing


def update_peptide_pool(seq:aa.AminoAcidSeqRecord,
        peptide_pool:Set[aa.AminoAcidSeqDict],
        label_counter:Set[aa.AminoAcidSeqRecord], label:str,
        update_label:bool=True) -> None:
    """ Add a peptide sequence to a peptide pool if not already exists,
    otherwise update the same peptide's description """
    if 'X' in seq.seq:
        return
    if update_label:
        if label not in label_counter:
            label_counter[label] = 0
        label_counter[label] += 1
        seq.description = label + '|' + str(label_counter[label])
        seq.id = seq.description
        seq.name = seq.description
        same_peptide = get_equivalent(peptide_pool, seq)
        if same_peptide:
            same_peptide:aa.AminoAcidSeqRecord
            same_peptide.description += ' ' + seq.description
            same_peptide.id = same_peptide.description
            same_peptide.name = same_peptide.description
        else:
            peptide_pool.add(seq)
    else:
        if label not in label_counter:
            label_counter[label] = 0
        label_counter[label] += 1
        seq.description = label + '|'\
            + str(label_counter[label])
        peptide_pool.add(seq)

class PVGPeptideFinder():
    """ Variant peptide pool as dict.

    Attributes:
        - `tx_id` (str): Transcript ID.
        - `peptides` (Dict[AminoAcidSeqRecord, Set[VariantPeptideMetadata]]):
          The peptide data pool, with keys being the AminoAcidRecord, and
          values being set of VariantPeptdieMetadata (variants and ORF).
        - `labels` (Dict[str,int]): Label counter, as a dict with key being the
          variant peptide label, and values being the number of occurrence
          of this label.
        - `global_variant` (VariantRecord): A variant that should be applied to
          every node. This is useful for PVG derived from a circRNA, and this
          argument is the circRNA variant.
    """
    PEPTIDE_FINDING_MODES = ['misc', 'archipel']

    def __init__(self, tx_id:str, peptides:TypeVariantPeptideMetadataMap=None,
            seqs:Set[Seq]=None, labels:Dict[str,int]=None, mode:str='misc',
            global_variant:VariantRecord=None, gene_id:str=None,
            truncate_sec:bool=False, w2f:bool=False, check_external_variants:bool=True,
            cleavage_params:CleavageParams=None, check_orf:bool=False):
        """ constructor """
        assert mode in self.PEPTIDE_FINDING_MODES
        self.tx_id = tx_id
        self.peptides = peptides or {}
        self.mode = mode
        self.seqs = seqs or set()
        self.labels = labels or {}
        self.global_variant = global_variant
        self.gene_id = gene_id
        self.truncate_sec = truncate_sec
        self.w2f = w2f
        self.check_external_variants = check_external_variants
        self.cleavage_params = cleavage_params
        self.check_orf = check_orf

    def find_candidate_node_paths_misc(self, node:PVGNode, orfs:List[PVGOrf],
            cleavage_params:CleavageParams, tx_id:str, gene_id:str,
            leading_node:PVGNode, subgraphs:SubgraphTree, is_circ_rna:bool,
            backsplicing_only:bool
            ) -> PVGCandidateNodePaths:
        """ Find all miscleaved nodes.

        node vs leading_node:
            - When finding candidate peptides sthat starts with a series of
              contigiuos nodes, a copy of the first node must be first created.
              This is because if the ORF start site is located in the node, the
              valid peptide sequence should be called form the ORF start site,
              rather than the beginning of the sequence. So then a copy should
              be created from this node and only keep the sequence after the ORF start.
            - The `leading_node` must tbe the original copy of the first node
              of the node path in the graph, that the first element of the `node_path`
              is copied from. This is used to retrieve INDELs from the downstream
              node.

        Args:
            - `node` (PVGNode): The node that miscleaved peptide sequences
              starts from it should be called.
            - `orfs` (List[PVGPOrf]): The ORF start and end locations.
            - `cleavage_params` (CleavageParams): Cleavage related parameters.
            - `tx_id` (str): Transcript ID.
            - `gene_id` (str): Gene ID.
            - `leading_nodes` (List[PVGNode]): The leading node to find candidate
              node paths.
            - `leading_node` (PVGNode): The start node that the miscleaved
              peptides are called from. This node must present in the PVG graph.
        """
        if not orfs:
            raise ValueError('ORFs are empty')
        cur_path = PVGNodePath([node], set())
        queue = deque([cur_path])
        paths = PVGCandidateNodePaths(
            data=deque([]),
            cleavage_params=cleavage_params,
            orfs=orfs,
            tx_id=tx_id,
            gene_id=gene_id,
            leading_node=leading_node,
            subgraphs=subgraphs,
            is_circ_rna=is_circ_rna
        )

        if not (node.cpop_collapsed or node.truncated) and \
                (not backsplicing_only or len(node.get_subgraph_id_set()) > 1):
            additional_variants = leading_node.get_downstream_stop_altering_variants()
            path = PVGNodePath([node], additional_variants)
            if path.is_too_long(self.cleavage_params) and not node.selenocysteines:
                return paths
            if not path.is_too_short(self.cleavage_params):
                paths.data.append(path)

        while queue:
            cur_path = queue.pop()
            cur_node = cur_path.nodes[-1]
            # Turn it into a doct of id to variant
            path_vars = cur_path.variants

            n_cleavages = len([x for x in cur_path.nodes if not x.cpop_collapsed]) - 1

            if n_cleavages >= cleavage_params.miscleavage:
                continue

            # This is done to reduce the complexity when the node has too
            # many out nodes, and its out nodes also have too many out nodes.
            if cleavage_params.additional_variants_per_misc == -1:
                allowed_n_vars = float('Inf')
            else:
                allowed_n_vars = cleavage_params.max_variants_per_node
                allowed_n_vars += (n_cleavages + 1) * cleavage_params.additional_variants_per_misc

            for _node in cur_node.out_nodes:
                if is_circ_rna and _node.is_hybrid_node(subgraphs):
                    continue
                if _node.truncated:
                    continue
                new_path = copy.copy(cur_path)

                is_stop = len(_node.seq.seq) == 1 and _node.seq.seq.startswith('*')
                if is_stop:
                    continue

                new_path.append(_node)

                if backsplicing_only and not cur_path.is_subgraph_spanning():
                    # Skip peptides not spaning the backsplicing site
                    continue

                additional_variants = _node.get_downstream_stop_altering_variants()

                cur_vars = set(path_vars.keys())
                for var in _node.variants:
                    cur_vars.add(var.variant.id)
                cur_vars.update({v.id for v in additional_variants})
                if len(cur_vars) > allowed_n_vars:
                    continue

                if not _node.cpop_collapsed:
                    new_path.add_additional_variants(additional_variants)
                    if new_path.is_too_long(self.cleavage_params) \
                            and not new_path.has_trailing_selenocysteins():
                        continue
                    if not new_path.is_too_short(self.cleavage_params):
                        paths.data.append(new_path)
                    if n_cleavages + 1 == cleavage_params.miscleavage:
                        continue
                    if n_cleavages + 1 > cleavage_params.miscleavage:
                        raise ValueError('Something just went wrong')
                queue.append(new_path)
        return paths

    def find_candidate_node_paths_archipel(self, node:PVGNode,
            orfs:List[PVGOrf], cleavage_params:CleavageParams, tx_id:str, gene_id:str,
            leading_node:PVGNode, subgraphs:SubgraphTree, is_circ_rna:bool,
            backsplicing_only:bool, is_start_codon:bool):
        """
        This function is used to find node paths in two styles: archipel and reef.
        Archipel nodes are varaint nodes as islands, separated by reference nodes,
        and are surrounded by flanking reference nodes (e.g., the ocean). Reef
        nodes are those only carrying backbone alterning variants, such as large
        indels (alternative splicing), fusion, and circRNA, without any other
        small variants.
        """
        if not orfs:
            raise ValueError('ORFs are empty')
        if any(v.variant != node.global_variant for v in node.variants):
            islands = [0]
            nflanking = []
        else:
            islands = []
            nflanking = [0]
        flanking_size = cleavage_params.flanking_size
        cur_path = PVGNodePathArchipel(
            nodes=[node], additional_variants=set(),
            nflanking=nflanking, flanking_size=flanking_size,
            islands=islands
        )

        queue:Deque[PVGNodePath] = deque([cur_path])
        paths = PVGCandidateNodePaths(
            data=deque([]),
            cleavage_params=cleavage_params,
            orfs=orfs,
            tx_id=tx_id,
            gene_id=gene_id,
            leading_node=leading_node,
            subgraphs=subgraphs,
            is_circ_rna=is_circ_rna
        )

        while queue:
            cur_path = queue.pop()

            for out_node in cur_path.nodes[-1].out_nodes:
                is_stop = len(out_node.seq.seq) == 1 and out_node.seq.seq.startswith('*')
                if is_stop:
                    if cur_path.is_valid_path() \
                            and (not backsplicing_only or new_path.is_subgraph_spanning()):
                        paths.data.append(cur_path)
                    continue

                if isinstance(cur_path, PVGNodePathArchipel):
                    new_path = copy.copy(cur_path)
                    i = new_path.number_of_nodes()
                    new_path.append(out_node)
                    if not cur_path.has_any_potential_island_variant(out_node):
                        if not new_path.has_any_island():
                            if new_path.is_nflanking_full():
                                if out_node.global_variant:
                                    new_path = new_path.to_reef()
                                    queue.append(new_path)
                                continue
                            new_path.nflanking.append(i)
                        else:
                            if len(new_path.cflanking) > flanking_size:
                                continue
                            new_path.cflanking.append(i)
                            if len(new_path.cflanking) == flanking_size:
                                if not backsplicing_only or new_path.is_subgraph_spanning():
                                    paths.data.append(new_path)
                                continue
                    else:
                        if not is_start_codon and len(new_path.nflanking) < flanking_size:
                            continue
                        new_path.islands.append(i)
                        new_path.cflanking = []
                elif isinstance(cur_path, PVGNodePathReef):
                    if cur_path.has_any_potential_island_variant(out_node):
                        paths.data.append(cur_path)
                        continue
                    if cur_path.is_last_reef_full():
                        paths.data.append(cur_path)
                        continue
                    new_path = copy.copy(cur_path)
                    new_path.add_node(out_node)
                else:
                    raise TypeError(f"Unrecognized type {type(cur_path)}")
                queue.append(new_path)
        return paths

    def add_peptide_sequences(self, node:PVGNode, orfs:List[PVGOrf],
            cleavage_params:CleavageParams, check_variants:bool, is_start_codon:bool,
            additional_variants:List[VariantRecord], denylist:Set[str],
            leading_node:PVGNode=None, subgraphs:SubgraphTree=None,
            circ_rna:CircRNAModel=None, backsplicing_only:bool=False):
        """ Add amino acid sequences starting from the given node, with number
        of miscleavages no more than a given number. The sequences being added
        are the sequence of the current node, and plus n of downstream nodes,
        with n ranges from 1 to miscleavages.

        Args:
        - `node` (PVGNode): The start node.
        - `orfs` (List[int, int]): The start and end position of the ORF.
        - `cleavage_params` (CleavageParams): Cleavage related parameters.
        - `check_variants` (bool): Whether to check for variants.
        - `is_start_codon` (bool): Whether the peptide starts with the start
            codon (M).
        - `additional_variants` (List[VariantRecord]): Additional variants,
            e.g., start gain, cleavage gain, stop lost.
        - `denylist` (Set[str]): Peptide sequences that should be excluded.
        - `leading_node` (PVGNode): The start node that the miscleaved
            peptides are called from. This node must be present in the PVG graph.
        - `subgraphs` (SubgraphTree): The subgraph tree that holds how subgraphs
            are incorporated into the main graph.
        - `circ_rna` (CircRNAModel): The circRNA molecular model from which variant
            peptides are being called.
        - `backsplicing_only` (bool): Whether to only output variant peptides
            spanning the backsplicing site.
        """
        if leading_node is None:
            leading_node = node

        if self.mode == 'misc':
            candidate_paths = self.find_candidate_node_paths_misc(
                node=node, orfs=orfs, cleavage_params=cleavage_params,
                tx_id=self.tx_id, gene_id=self.gene_id, leading_node=leading_node,
                subgraphs=subgraphs, is_circ_rna=circ_rna is not None,
                backsplicing_only=backsplicing_only
            )
        else:
            candidate_paths = self.find_candidate_node_paths_archipel(
                node=node, orfs=orfs, cleavage_params=cleavage_params,
                tx_id=self.tx_id, gene_id=self.gene_id, leading_node=leading_node,
                subgraphs=subgraphs, is_circ_rna=circ_rna is not None,
                backsplicing_only=backsplicing_only, is_start_codon=is_start_codon
            )
        if self.global_variant and self.global_variant not in additional_variants:
            additional_variants.append(self.global_variant)

        it = candidate_paths.join_peptides(
            pool=self.peptides,
            check_variants=check_variants,
            additional_variants=additional_variants,
            denylist=denylist,
            is_start_codon=is_start_codon,
            circ_rna=circ_rna,
            truncate_sec=self.truncate_sec,
            check_external_variants=self.check_external_variants,
            check_orf=self.check_orf
        )
        for seq, metadata in it:
            if 'X' in seq:
                continue
            if '*' in seq:
                raise ValueError('Invalid amino acid symbol found in the sequence.')
            val = self.peptides.setdefault(seq, {})
            key = metadata.get_key()
            if key not in val:
                val[key] = metadata
            self.seqs.add(seq)

    def find_codon_reassignments(self, seq:Seq, w2f:bool=False) -> List[VariantRecord]:
        """ Find potential codon reassignments from the given sequence. """
        variants = []
        if w2f:
            i = 0
            while i > -1:
                i = seq.find('W', start=i)
                if i > -1:
                    w2f_var = create_variant_w2f(self.tx_id, i)
                    variants.append(w2f_var)
                    i += 1
                    if i >= len(seq):
                        break
        return variants

    def is_valid_seq(self, seq:Seq, denylist:Set[str]) -> bool:
        """ Check whether the seq is valid """
        if seq in self.seqs:
            return True
        min_mw = self.cleavage_params.min_mw
        return self.seq_has_valid_size(seq) \
            and seq not in denylist \
            and 'X' not in seq \
            and SeqUtils.molecular_weight(seq, 'protein') >= min_mw

    def seq_has_valid_size(self, seq:Seq) -> bool:
        """ Check whether the seq has valid length """
        min_length = self.cleavage_params.min_length
        max_length = self.cleavage_params.max_length

        return min_length <= len(seq) <= max_length

    def translational_modification(self, w2f:bool, denylist:Set[str]):
        """ Sequence level translational modification, e.g., W2F reassignment """
        for seq in copy.copy(self.peptides):
            reassignments:List[seqvar.VariantRecord] = []
            reassignments += self.find_codon_reassignments(seq, w2f)
            if not reassignments:
                continue

            # W > F codon reassignments
            for k in range(1, len(reassignments) + 1):
                for comb in itertools.combinations(reassignments, k):
                    seq_mod = copy.copy(seq)
                    for v in comb:
                        seq_mod_new = seq_mod[:v.location.start] + v.alt
                        if v.location.end < len(seq):
                            seq_mod_new += seq_mod[v.location.end:]
                        seq_mod = seq_mod_new

                    if not self.is_valid_seq(seq_mod, denylist):
                        continue
                    for metadata in self.peptides[seq].values():
                        cur_metadata = copy.copy(metadata)
                        cur_metadata.label += '|' + '|'.join(v.id for v in comb)
                        cur_metadata.has_variants = True

                        val = self.peptides.setdefault(seq_mod, {})
                        key = cur_metadata.get_key()
                        if key not in val:
                            val[key] = cur_metadata
                        self.seqs.add(seq_mod)

    def get_peptide_sequences(self, keep_all_occurrence:bool=True,
            orf_id_map:Dict[Tuple[int,int],str]=None, check_variants:bool=True
            ) -> Dict[Seq,List[AnnotatedPeptideLabel]]:
        """ Get the peptide sequence and check for miscleavages.

        Args:
            - `keep_all_occurrence` (bool): Whether to keep all occurance of
              the peptide within the graph/transcript.
            - `orf_id_map` (Dict[int, str]): An ID map of ORF IDs from ORF
              start position.
        """
        peptide_segments: Dict[Seq, List[AnnotatedPeptideLabel]] = {}

        for seq in self.peptides:
            if '*' in seq:
                raise ValueError('Invalid amino acid symbol found in the sequence.')
            labels = []
            metadatas = list(self.peptides[seq].values())
            metadatas.sort(key=lambda x: x.orf[0])
            unique_labels = set()
            had_pure_circ_rna = False
            for metadata in metadatas:
                if check_variants and not metadata.has_variants:
                    continue
                if metadata.is_pure_circ_rna:
                    if had_pure_circ_rna:
                        continue
                    had_pure_circ_rna = True
                orf_id = None
                if orf_id_map:
                    orf_id = orf_id_map[metadata.orf]

                label = metadata.label
                if orf_id:
                    label += f"|{orf_id}"

                if label in unique_labels:
                    continue
                unique_labels.add(label)

                if label in self.labels:
                    self.labels[label] += 1
                else:
                    self.labels[label] = 1
                label += f"|{self.labels[label]}"
                labels.append(label)
                seq_anno = AnnotatedPeptideLabel(label, metadata.segments)
                if seq in peptide_segments:
                    peptide_segments[seq].append(seq_anno)
                else:
                    peptide_segments[seq] = [seq_anno]

                if not keep_all_occurrence:
                    break

            if not labels:
                continue

        return peptide_segments
