""" SubgraphMapper """
from __future__ import annotations
from typing import Dict, Tuple, TYPE_CHECKING
import uuid
from moPepGen.SeqFeature import FeatureLocation


if TYPE_CHECKING:
    from moPepGen.seqvar import VariantRecord

class SubgraphBranch():
    """ ParentGraphLocation """
    def __init__(self, _id:str, level:int, parent_id:str,  location:FeatureLocation,
            variant:VariantRecord, feature_type:str, feature_id:str):
        """ constructur """
        self.id = _id
        self.level = level
        self.parent_id = parent_id
        self.location = location
        self.variant = variant
        self.feature_type = feature_type
        self.feature_id= feature_id

    def is_parent(self, other:SubgraphBranch) -> bool:
        """ Checks if it is the parent of a given subgraph """
        return self.id == other.parent_id

class SubgraphTree():
    """ This SubgraphTree class defines the relationship of subgraphes and
    parents/children in a tree-like data structure. """
    def __init__(self, data:Dict[str, SubgraphBranch]=None,
            root:SubgraphBranch=None):
        """ Constructor """
        self.data = data or {}
        self.root = root

    def __getitem__(self, key:str) -> SubgraphBranch:
        """ Get item """
        return self.data[key]

    def add_root(self, _id:str, feature_type:str, feature_id:str, variant:VariantRecord):
        """ Add the main graph """
        subgraph = SubgraphBranch(
            _id=_id, level=0, parent_id=None, location=None, variant=variant,
            feature_type=feature_type, feature_id=feature_id
        )
        self.data[_id] = subgraph
        self.root = subgraph

    def add_subgraph(self, child_id:str, parent_id:str, level:int,
            start:int, end:int, variant:VariantRecord,
            feature_type:str, feature_id:str):
        """ Add subgraph """
        location = FeatureLocation(
            seqname=parent_id, start=start, end=end
        )
        subgraph = SubgraphBranch(
            _id=child_id, level=level, parent_id=parent_id,
            location=location, variant=variant,
            feature_type=feature_type, feature_id=feature_id
        )
        self.data[child_id] = subgraph

    @staticmethod
    def generate_subgraph_id() -> str:
        """ Generate a random subgraph ID """
        return str(uuid.uuid4())

    def find_compatible_parents(self, first:str, second:str
            ) -> Tuple[SubgraphBranch, SubgraphBranch]:
        """ For two given subgraphs, find their parent graphs that are
        compatible. Compatible means either they share the same parent subgraph,
        or one is the parent of the other subgraph. """
        subgraph1 = self[first.subgraph_id]
        subgraph2 = self[second.subgraph_id]

        while subgraph1.id != subgraph2.id and \
                subgraph1.parent_id != subgraph2.parent_id and \
                subgraph1.id !=  subgraph2.parent_id and \
                subgraph1.parent_id != subgraph2.id:
            level1 = subgraph1.level
            level2 = subgraph2.level
            if level1 >= level2:
                subgraph1 = self[subgraph1.parent_id]
            if level2 >= level1:
                subgraph2 = self[subgraph2.parent_id]

        return subgraph1, subgraph2
