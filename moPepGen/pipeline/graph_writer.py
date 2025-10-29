"""Graph serialization utilities for the variant peptide calling pipeline.

This module provides the GraphWriter class for saving variant graphs (DNA and
protein level) to JSON files. This is useful for debugging, visualization, and
quality control of the graph construction process.

The graphs are saved with structured filenames that identify the transcript,
variant source (main/fusion/circRNA), and graph type (TVG/CVG/PVG).
"""
from __future__ import annotations
from typing import Dict
import json
from pathlib import Path
from .models import DGraphs, PGraphs

class GraphWriter:
    """Handles serialization of variant graphs to JSON files.

    This class provides methods to save DNA-level graphs (TVG/CVG) and
    protein-level graphs (PVG) for different variant sources. The graphs
    are saved with structured filenames that enable easy identification
    and retrieval.

    The file naming convention is:
        {tx_id}_main_{graph_type}.json              # Main transcriptional variants
        {tx_id}_Fusion_{variant_id}_{graph_type}.json   # Fusion variants
        {tx_id}_circRNA_{variant_id}_{graph_type}.json  # CircRNA variants

    Where graph_type is one of: TVG (Three-frame Transcript Variant Graph),
    CVG (Chimeric Variant Graph for fusions/circRNA), or PVG (Peptide Variant Graph).

    Attributes:
        out_dir: Directory where graph JSON files will be written. Created
            automatically if it doesn't exist.
    """
    def __init__(self, out_dir: Path):
        """Initialize the GraphWriter.

        Args:
            out_dir: Output directory for graph JSON files. Will be created
                if it doesn't exist (including parent directories).
        """
        self.out_dir = out_dir

    def _write_json(self, file_path: Path, data: Dict):
        """Write a dictionary to a JSON file.

        Creates parent directories if they don't exist. The JSON is written
        in a human-readable format (indented).

        Args:
            file_path: Full path where the JSON file should be written
            data: Dictionary to serialize as JSON. Typically the result of
                calling the graph's jsonfy() method.
        """
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, 'wt') as handle:
            json.dump(data, handle)

    def write_dgraphs(self, tx_id: str, dgraphs: DGraphs):
        """Write DNA-level graphs (TVG/CVG) to JSON files.

        Saves the main transcript variant graph and any fusion or circRNA
        graphs to separate JSON files with structured names.

        Args:
            tx_id: Transcript identifier (e.g., 'ENST00000123456.2')
            dgraphs: Tuple of DNA-level graphs in the format:
                - dgraphs[0]: Main TVG (ThreeFrameTVG) or None
                - dgraphs[1]: Dict mapping fusion IDs to their CVG graphs
                - dgraphs[2]: Dict mapping circRNA IDs to their CVG graphs

        Output files:
            - {tx_id}_main_TVG.json: Main transcriptional variants
            - {tx_id}_Fusion_{variant_id}_TVG.json: Per-fusion graph
            - {tx_id}_circRNA_{variant_id}_TVG.json: Per-circRNA graph

        Note:
            Only non-None graphs are written. Empty graphs are skipped.
        """
        if dgraphs[0]:
            self._write_json(self.out_dir / f"{tx_id}_main_TVG.json", dgraphs[0].jsonfy())
        for var_id, graph in dgraphs[1].items():
            if graph is None:
                continue
            self._write_json(self.out_dir / f"{tx_id}_Fusion_{var_id}_TVG.json", graph.jsonfy())
        for var_id, graph in dgraphs[2].items():
            if graph is None:
                continue
            self._write_json(self.out_dir / f"{tx_id}_circRNA_{var_id}_TVG.json", graph.jsonfy())

    def write_pgraphs(self, tx_id: str, pgraphs: PGraphs):
        """Write protein-level graphs (PVG) to JSON files.

        Saves the main peptide variant graph and any fusion or circRNA
        peptide graphs to separate JSON files with structured names.

        Args:
            tx_id: Transcript identifier (e.g., 'ENST00000123456.2')
            pgraphs: Tuple of protein-level graphs in the format:
                - pgraphs[0]: Main PVG (PeptideVariantGraph) or None
                - pgraphs[1]: Dict mapping fusion IDs to their PVGs
                - pgraphs[2]: Dict mapping circRNA IDs to their PVGs

        Output files:
            - {tx_id}_main_PVG.json: Main peptide variant graph
            - {tx_id}_Fusion_{variant_id}_PVG.json: Per-fusion peptide graph
            - {tx_id}_circRNA_{variant_id}_PVG.json: Per-circRNA peptide graph

        Note:
            Only non-None graphs are written. Empty graphs are skipped.
            These graphs represent the translated (protein-level) versions
            of the DNA graphs, after applying translation and cleavage rules.
        """
        if pgraphs[0]:
            self._write_json(self.out_dir / f"{tx_id}_main_PVG.json", pgraphs[0].jsonfy())
        for var_id, graph in pgraphs[1].items():
            if graph is None:
                continue
            self._write_json(self.out_dir / f"{tx_id}_Fusion_{var_id}_PVG.json", graph.jsonfy())
        for var_id, graph in pgraphs[2].items():
            if graph is None:
                continue
            self._write_json(self.out_dir / f"{tx_id}_circRNA_{var_id}_PVG.json", graph.jsonfy())
