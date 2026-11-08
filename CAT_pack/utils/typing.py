"""Shared type definitions for CAT_pack."""

from enum import StrEnum
from typing import Hashable, TypeVar

Taxid = TypeVar("Taxid", bound=Hashable)
"""Type of a taxonomic identifier in a lineage"""

class EntityType(StrEnum):
    CONTIG = "contig"
    BIN = "bin"
    SEQUENCE = "sequence"
