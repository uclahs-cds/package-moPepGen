from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple, Any

# We keep types loosely annotated to avoid heavy imports at module import time.

class GraphPhase(str, Enum):
    TVG = "TVG"
    PVG = "PVG"

@dataclass
class TimeoutTracer:
    backbone: Optional[str] = None  # 'main' | 'fusion' | 'circRNA'
    graph: Optional[GraphPhase] = None  # GraphPhase

@dataclass
class Flags:
    noncanonical_transcripts: bool
    backsplicing_only: bool
    coding_novel_orf: bool
    skip_failed: bool
    truncate_sec: bool
    w2f_reassignment: bool

@dataclass
class Limits:
    max_adjacent_as_mnv: int
    max_variants_per_node: List[int] = field(default_factory=lambda: [7])
    additional_variants_per_misc: List[int] = field(default_factory=lambda: [2])
    in_bubble_cap_step_down: int = 0

# Result graph container types: (main, {fusion_id: graph}, {circ_id: graph})
DGraphs = Tuple[Any, Dict[str, Any], Dict[str, Any]]
PGraphs = Tuple[Any, Dict[str, Any], Dict[str, Any]]

@dataclass
class CallVariantDispatch:
    tx_id: str
    variant_series: Any
    tx_seqs: Dict[str, Any]
    gene_seqs: Dict[str, Any]
    reference_data: Any
    codon_table: Any
    pool: Any
    cleavage_params: Any
    flags: Flags
    limits: Limits
    save_graph: bool
    timeout_seconds: int

@dataclass
class CallResult:
    tx_id: str
    peptide_anno: Dict[Any, List[Any]]  # Dict[Seq, List[AnnotatedPeptideLabel]]
    dgraphs: DGraphs
    pgraphs: PGraphs
    success: Dict[str, bool]  # keys: 'variant', 'fusion', 'circRNA'
