""" Module for variant peptide dict """
from __future__ import annotations
from collections import deque
import copy
from typing import Deque, Dict, Iterable, List, Set, Tuple, TYPE_CHECKING
from moPepGen import aa, get_equivalent
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.aa import VariantPeptideIdentifier as vpi


if TYPE_CHECKING:
    from moPepGen.seqvar.VariantRecord import VariantRecord

class MiscleavedNodes():
    """ Helper class for looking for peptides with miscleavages """
    def __init__(self, data:Deque[List[PVGNode]], orf:Tuple[int,int]=None):
        """ constructor """
        self.data = data
        self.orf = orf

    @staticmethod
    def find_miscleaved_nodes(node:PVGNode, orf:Tuple[int,int], miscleavage:int
            ) -> MiscleavedNodes:
        """ find all miscleaved nodes """
        if orf[0] is None:
            raise ValueError('orf is None')
        queue = deque([[node]])
        nodes = MiscleavedNodes(deque([[node]]), orf)
        i = 0
        while i < miscleavage and queue:
            i += 1
            next_cleavage = deque([])
            while queue:
                cur_batch = queue.pop()
                cur_node = cur_batch[-1]
                for _node in cur_node.out_nodes:
                    if not _node.out_nodes:
                        continue
                    if _node.truncated:
                        continue
                    new_batch = copy.copy(cur_batch)
                    stop_index = _node.seq.seq.find('*')
                    if stop_index > -1:
                        left_node = _node[:stop_index + 1]
                        left_node.seq = left_node.seq[:stop_index]
                        new_batch.append(left_node)
                        nodes.data.append(new_batch)
                        continue
                    new_batch.append(_node)
                    next_cleavage.appendleft(new_batch)
                    nodes.data.append(copy.copy(new_batch))
            queue = next_cleavage
        return nodes

    def join_miscleaved_peptides(self, check_variants:bool,
            additional_variants:List[VariantRecord], is_start_codon:bool=False
            ) -> Iterable[Tuple[aa.AminoAcidSeqRecord, VariantPeptideMetadata]]:
        """ join miscleaved peptides and update the peptide pool.

        Args:
            peptide_pool (Set[aa.AminoAcidSeqRecord]): The container for all
                peptides called from the PeptdieVariantGraph.
            graph (PeptideVariantGraph): The graph object.
            check_variants (bool): When true, only peptides that carries at
                least 1 variant are kept. And when false, all unique peptides
                are reported （e.g. noncoding).
            label_counter (Dict[str,int]): The object counts the total
                occurrences of each variant label. An int number is appended
                to the end of the label (e.g. ENST0001|SNV-10-T-C|5)
        """
        for queue in self.data:
            metadata = VariantPeptideMetadata(orf=self.orf)
            seq:aa.AminoAcidSeqRecord = None

            for node in queue:
                other = copy.deepcopy(node.seq)
                other.__class__ = aa.AminoAcidSeqRecord
                other:aa.AminoAcidSeqRecord
                if seq is None:
                    seq = other
                else:
                    seq = seq + other
                if check_variants:
                    # metadata.variants.update(node.frameshifts)
                    for variant in node.variants:
                        metadata.variants.add(variant.variant)
                    metadata.variants.update(additional_variants)

            if seq:
                seq.__class__ = aa.AminoAcidSeqRecord

            if check_variants and not metadata.variants:
                continue

            yield seq, metadata

            if is_start_codon and seq.seq.startswith('M'):
                seq = seq[1:]
                yield seq, metadata

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

class VariantPeptideMetadata():
    """ Variant peptide metadata """
    def __init__(self, variants:Set[VariantRecord]=None,
            orf:Tuple[int,int]=None):
        """  """
        self.variants = variants or set()
        self.orf = orf

T = Dict[aa.AminoAcidSeqRecord, Set[VariantPeptideMetadata]]

class VariantPeptideDict():
    """ Variant peptide pool as dict.

    Attributes:
        tx_id (str): Transcript ID.
        peptides (Dict[aa.AminoAcidSeqRecord, Set[VariantPeptideMetadata]]):
            The peptide data pool, with keys being the AminoAcidRecord, and
            values being set of VariantPeptdieMetadata (variants and ORF).
        labels (Dict[str,int]): Label counter, as a dict with key being the
            variant piptide label, and values being the number of occurrence
            of this lable.
    """
    def __init__(self, tx_id:str, peptides:T=None, labels:Dict[str,int]=None):
        """ constructor """
        self.tx_id = tx_id
        self.peptides = peptides or {}
        self.labels = labels or {}

    def add_miscleaved_sequences(self, node:PVGNode, orf:List[int, int],
            miscleavage:int, check_variants:bool, is_start_codon:bool,
            additional_variants:List[VariantRecord]):
        """ Add amino acid sequences starting from the given node, with number
        of miscleavages no more than a given number. The sequences being added
        are the sequence of the current node, and plus n of downstream nodes,
        with n ranges from 1 to miscleavages.

        Args:
            node (PVGNode): The start node.
            miscleavage (int): Number of miscleavage to allow.
            check_variants (bool): Whether to check variants.
            has_known_orf (bool): Whether has known ORF. False for noncoding
                genes/transcripts.
        """
        miscleaved_nodes = MiscleavedNodes.find_miscleaved_nodes(
            node=node, orf=orf, miscleavage=miscleavage
        )
        for seq, metadata in miscleaved_nodes.join_miscleaved_peptides(
                check_variants, additional_variants, is_start_codon):
            if 'X' in seq.seq:
                continue
            if '*' in seq:
                raise ValueError('Invalid amino acid symbol found in the sequence.')
            val = self.peptides.setdefault(seq, set())
            val.add(metadata)

    def get_peptide_sequences(self, keep_all_occurrence:bool=True,
            orf_id_map:Dict[Tuple[int,int],str]=None
            ) -> Set[aa.AminoAcidSeqRecord]:
        """ Get the peptide sequence and check for miscleavages.

        Args:
            keep_all_occurrence (bool): Whether to keep all occurance of the
                peptide within the graph/transcript.
            orf_id_map (Dict[int, str]): An ID map of ORF IDs from ORF start
                position.
        """
        peptide_pool:Set[aa.AminoAcidSeqRecord] = set()

        for seq, metadatas in self.peptides.items():
            if '*' in seq:
                raise ValueError('Invalid amino acid symbol found in the sequence.')
            labels = []
            metadatas = list(metadatas)
            metadatas.sort(key=lambda x: x.orf[0])
            has_pure_circ_rna = False
            for metadata in metadatas:
                variants = list(metadata.variants)
                orf_id = None
                if orf_id_map:
                    orf_id = orf_id_map[metadata.orf]

                label = vpi.create_variant_peptide_id(self.tx_id, variants, orf_id)

                is_pure_circ_rna = len(variants) == 1 and variants[0].is_circ_rna()
                if is_pure_circ_rna:
                    if has_pure_circ_rna:
                        continue
                    has_pure_circ_rna = True

                if label in self.labels:
                    self.labels[label] += 1
                else:
                    self.labels[label] = 1
                label += f"|{self.labels[label]}"
                labels.append(label)

                if not keep_all_occurrence:
                    break
            seq.description = " ".join(labels)
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

        return peptide_pool
