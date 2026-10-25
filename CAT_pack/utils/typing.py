"""Shared type definitions for CAT_pack."""

from __future__ import absolute_import

from enum import Enum
from typing import Hashable, TypeVar

Taxid = TypeVar("Taxid", bound=Hashable)
"""Type of a taxonomic identifier in a lineage"""

class EntityType(str, Enum):
    CONTIG = "contig"
    BIN = "bin"
    SEQUENCE = "sequence"

    def __str__(self) -> str:
        return self.value