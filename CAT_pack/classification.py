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



@dataclass(slots=True)
class ORFClassification:
    orf_id: str
    status: ORFStatus

    n_hits: int = 0
    taxid: str | None = None # TODO: replace with int down the line
    top_bitscore: Decimal | None = None
    lineage: tuple[str, ...] = () # TODO: update on lineage internal change

