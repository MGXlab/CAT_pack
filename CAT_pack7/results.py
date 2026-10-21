from dataclasses import dataclass
from enum import auto, IntEnum

import numpy as np

from .io.parsers import AlignmentInput
from .tax import Taxonomy, WeightedLCA


@dataclass(frozen=True, slots=True)
class Lineage:
    taxids: tuple[str, ...]

    def __str__(self) -> str:
        return ";".join(reversed(self.taxids))

    def format(self, branches: set[str], no_stars: bool = False) -> str:
        if no_stars or len(self.taxids) <= 2:
            return str(self)
        taxids = list(self.taxids)
        for index, parent in enumerate(self.taxids[1:]):
            if parent in branches:
                break
            taxids[index] += "*"
        return ";".join(reversed(taxids))


class ClassificationStatus(IntEnum):
    ASSIGNED = auto()
    NO_ORFS = auto()
    NO_HITS = auto()
    NO_TAXIDS = auto()
    NO_LINEAGE_SUPPORT = auto()

    @property
    def reason(self) -> str:
        return {
            ClassificationStatus.NO_ORFS: "no ORFs found",
            ClassificationStatus.NO_HITS: "no hits to database",
            ClassificationStatus.NO_TAXIDS: "hits not found in taxonomy files",
            ClassificationStatus.NO_LINEAGE_SUPPORT: "no lineage reached minimum bit-score support",
        }[self]


@dataclass(frozen=True)
class Classification:
    """Classification of all contigs or bins "entities"
    stored as 1d numpy arrays"""
    entity_type: str
    entities: list[str]
    status: np.ndarray
    orf_offsets: np.ndarray
    orf_names: list[str]
    orf_row: np.ndarray
    orf_lca: np.ndarray
    orf_n_hits: np.ndarray
    orf_top_bitscore: np.ndarray
    n_orfs_with_taxid: np.ndarray
    assignment_offsets: np.ndarray
    weighted: WeightedLCA
    alignment: AlignmentInput
    taxonomy: Taxonomy

    def __len__(self) -> int:
        return len(self.entities)

    def get_status(self, entity: int) -> ClassificationStatus:
        return ClassificationStatus(int(self.status[entity]))

    def get_orfs(self, entity: int) -> slice:
        return slice(self.orf_offsets[entity], self.orf_offsets[entity + 1])

    def get_assignments(self, entity: int) -> range:
        return range(self.assignment_offsets[entity], self.assignment_offsets[entity + 1])

    @property
    def n_assigned(self) -> int:
        return int(np.count_nonzero(self.status == ClassificationStatus.ASSIGNED.value))
