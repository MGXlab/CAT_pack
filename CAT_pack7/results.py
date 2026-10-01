from dataclasses import dataclass
from decimal import Decimal
from enum import Enum, auto


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


class ORFStatus(Enum):
    """Classification status of a single predicted ORF.

    ASSIGNED: has taxid and LCA assigned to ORF
    NO_HIT: ORF has no accepted homology hit
    NO_TAXID: ORF accepted homology hits but no usable taxid association
    """
    NO_HIT = auto()
    NO_TAXID = auto()
    ASSIGNED = auto()


class ClassificationStatus(Enum):
    ASSIGNED = auto()
    NO_ORFS = auto()
    NO_HITS = auto()
    NO_TAXIDS = auto()
    NO_LINEAGE_SUPPORT = auto()


@dataclass(slots=True)
class ORFClassification:
    orf_id: str
    status: ORFStatus
    n_hits: int = 0
    taxid: str | None = None  # TODO: replace with int down the line
    message: str | None = None
    top_bitscore: Decimal | None = None
    lineage: tuple[str, ...] = ()  # TODO: update on lineage internal change


@dataclass(frozen=True, slots=True)
class TaxonomicAssignment:
    """One taxonomic assignment supported for an entity.

    If f < 0.5 this can cause multiple classifications for one contig.
    """
    taxid: str
    lineage: tuple[str, ...]
    lineage_scores: tuple[Decimal, ...]


class TaxonomyNamespace(Enum):
    """Support for multiple taxonomy namespaces in the future."""
    NCBI = auto()


@dataclass(slots=True)
class ClassificationResult:
    """A unified classification result for CAT and BAT."""
    entity_id: str
    status: ClassificationStatus
    taxonomy_namespace: TaxonomyNamespace = TaxonomyNamespace.NCBI
    assignments: tuple[TaxonomicAssignment, ...] = ()
    total_n_ORFs: int = 0
    based_on_n_ORFs: int = 0
    orf_results: tuple[ORFClassification, ...] = ()

    @property
    def assigned(self) -> bool:
        return self.status == ClassificationStatus.ASSIGNED

    @property
    def reason(self) -> str:
        if self.assigned:
            return f"based on {self.based_on_n_ORFs}/{self.total_n_ORFs} ORFs"
        return {
            ClassificationStatus.NO_ORFS: "no ORFs found",
            ClassificationStatus.NO_HITS: "no hits to database",
            ClassificationStatus.NO_TAXIDS: "hits not found in taxonomy files",
            ClassificationStatus.NO_LINEAGE_SUPPORT: "no lineage reached minimum bit-score support",
        }[self.status]
