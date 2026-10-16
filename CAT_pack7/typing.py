"""Shared type definitions for CAT_pack7."""

from typing import Hashable, TypeVar


Taxid = TypeVar("Taxid", bound=Hashable)
"""Type of a taxonomic identifier in a lineage"""
