#!/usr/bin/env python3

"""
Shared taxonomic classification engine for CAT and BAT.

It makes the logic that is currently duplicated between CAT and BAT centrally Accesable.

"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from enum import Enum, auto


class ORFStatus(Enum):
    """Classification status of a single predicted ORF

    Attributes:
    ASSIGNED: has taxid  and LCA assigned to ORF
    NO_HIT: ORF has no accepted homology hit
    NO_TAXID: ORF accepted homology hits but no usable taxid assosiation
    """
    NO_HIT = auto()
    NO_TAXID = auto()
    ASSIGNED = auto()

class ClassificationStatus(Enum):
    ASSINGED = auto()
    NO_ORFS = auto()
    NO_HITS = auto()
    NO_TAXIDS = auto()
    NO_LINEAGE_SUPPORT = auto()


@dataclass(slots=True)
class ORFClassification:
    orf_id: str
    status: ORFStatus

    n_hits: int = 0
    taxid: str | None = None # TODO: replace with int down the line
    top_bitscore: Decimal | None = None
    lineage: tuple[str, ...] = () # TODO: update on lineage internal change

@dataclass(frozen=True, slots=True)
class TaxonomicAssignment:
    """
    One taxonomic assignment supported for an entity
    if f < 0.5 can cause multiple classifications in one contig.
    """
    taxid: int
    support: float


class TaxonomyNamespace:
    pass


@dataclass(slots=True)
class ClassificationResult:
    entity_id: int
    #entity_type: EntityType #TODO: for later addition (other types like viral)
    status: ClassificationStatus
    #sequence_type: SequenceType #TODO: For later addition
    taxonomy_namespace: TaxonomyNamespace #Support for more then one taxid system
    assignments: tuple[TaxonomicAssignment, ...]




