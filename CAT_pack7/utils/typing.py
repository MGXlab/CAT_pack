"""Shared type definitions for CAT_pack7."""

from __future__ import absolute_import

from typing import Hashable, TypeVar

Taxid = TypeVar("Taxid", bound=Hashable)
"""Type of a taxonomic identifier in a lineage"""
