#!/usr/bin/env python3

"""
Shared taxonomic classification engine for CAT and BAT.

It makes the logic that is currently duplicated between CAT and BAT centrally Accesable.

"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from enum import Enum, auto
from typing import Mapping

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

class ClassificationEngine:
    """
    Shared ORF and WG  classifier for CAT and BAT

    NOTE: Not planning on implementing algorythm changes/optimizations here yet.
    First changing architecture to Object oriented and after that is done and
    results are identical the optimization can be done.
    """
    def __init__(self,*,taxid2parent:Mapping[str,str], fastaid2taxid: Mapping[str, str], fraction:Decimal) -> None:
        self.taxid2parent = taxid2parent
        self.fastaid2taxid = fastaid2taxid
        self.fraction = fraction

    def classify_orf(self):
        pass

    def classify_group(self):
        pass

