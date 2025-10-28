from __future__ import annotations
from typing import Dict
import json
from pathlib import Path
from .models import DGraphs, PGraphs

class GraphWriter:
    def __init__(self, out_dir: Path):
        self.out_dir = out_dir

    def _write_json(self, file_path: Path, data: Dict):
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, 'wt') as handle:
            json.dump(data, handle)

    def write_dgraphs(self, tx_id: str, dgraphs: DGraphs):
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
