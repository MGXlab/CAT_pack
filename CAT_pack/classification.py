#!/usr/bin/env python3

"""
Shared taxonomic classification engine for CAT and BAT.

It makes the logic that is currently duplicated between CAT and BAT centrally Accesable.

"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from enum import Enum, auto
from typing import Mapping, Sequence
import tax

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
    taxid: str | None = None # TODO: replace with int down the line
    top_bitscore: Decimal | None = None
    lineage: tuple[str, ...] = () # TODO: update on lineage internal change

@dataclass(frozen=True, slots=True)
class TaxonomicAssignment:
    """
    One taxonomic assignment supported for an entity
    if f < 0.5 can cause multiple classifications in one contig.
    """
    taxid: str
    lineage: tuple[str, ...]
    lineage_scores: tuple


@dataclass(slots=True)
class TaxonomyNamespace:
    NCBI = auto()


@dataclass(slots=True)
class ClassificationResult:
    entity_id: str
    #entity_type: EntityType #TODO: for later addition (other types like viral)
    status: ClassificationStatus
    #sequence_type: SequenceType #TODO: For later addition
    taxonomy_namespace: TaxonomyNamespace = TaxonomyNamespace.NCBI #Support for more then one taxid system
    assignments: tuple[TaxonomicAssignment, ...] = ()
    total_n_ORFs: int = 0
    based_on_n_ORFs: int = 0
    orf_results: tuple[ORFClassification, ...] = ()

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

    def classify_orf(self, orf_id: str, hits: Sequence[tuple[str, Decimal]] | None):
        if not hits:
            return ORFClassification(
                orf_id=orf_id,
                status=ORFStatus.NO_HIT
            )
        taxid, top_bitscore = tax.find_LCA_for_ORF(
            hits, self.fastaid2taxid, self.taxid2parent)

        if taxid.startswith("no taxid found"):
            return ORFClassification(
                orf_id=orf_id,
                status=ORFStatus.NO_TAXID,
                n_hits=len(hits),
                top_bitscore=top_bitscore,
                taxid=taxid # TODO: this will give back no taxid found to taxid what will be turned into
                # a int down the line will need to introduce a message system into this class
            )

        lineage = tax.find_lineage(taxid, self.taxid2parent)

        # TODO implement lineage starring, implementation here is not usefull for testing and will need to be
        # removed due to taxid becoming int in the future within this method

        return ORFClassification(
            orf_id=orf_id,
            status=ORFStatus.ASSIGNED,
            n_hits=len(hits),
            taxid=taxid,
            top_bitscore=top_bitscore,
            lineage=tuple(lineage)
        )

    def classify_group(self, *, entity_id: str, orf_ids: Sequence[str], orf2hits: Mapping[str, Sequence[tuple[str, Decimal]]]):
        """
        Classificaion of one contig or bin from ORFs
        """
        if not orf_ids: return ClassificationResult(
            entity_id=entity_id,
            status=ClassificationStatus.NO_ORFS,
        )
        orf_results= []
        lca_ORFs = []

        for orf_id in orf_ids:
            # get a ORF classification for one ORF and append it
            result = self.classify_orf(orf_id, orf2hits.get(orf_id))
            orf_results.append(result)

            if result.status == ORFStatus.NO_HIT:
                continue

            if result.status == ORFStatus.NO_TAXID:
                lca_ORFs.append((result.taxid, result.top_bitscore)) # TODO: add message to ORFclassification instead of misusing taxid
                continue

            lca_ORFs.append((result.taxid, result.top_bitscore))

        if not lca_ORFs:
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_HITS,
                assignments=(),
                total_n_ORFs=len(orf_ids),
                based_on_n_ORFs=0,
                orf_results=tuple(orf_results)
            )

        lineages, lineages_scores, based_on_n_orfs = tax.find_weighted_LCA(
            lca_ORFs,
            self.taxid2parent,
            self.fraction
        )

        if lineages == "no ORFs with taxids found.":
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_TAXIDS,
                assignments=(),
                total_n_ORFs=len(orf_ids),
                based_on_n_ORFs=0,
                orf_results=tuple(orf_results)
            )

        if lineages == "no lineage whitelisted.":
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_LINEAGE_SUPPORT,
                assignments=(),
                total_n_ORFs=len(orf_ids),
                based_on_n_ORFs=0,
                orf_results=tuple(orf_results)
            )

        assignments: list [TaxonomicAssignment] = []

        assignments = []

        for i, lineage in enumerate(lineages):
            assignments.append(
                TaxonomicAssignment(
                    taxid=lineage[0],
                    lineage=tuple(lineage),
                    lineage_scores=tuple(lineages_scores[i]),
                )
            )
        return ClassificationResult(
            entity_id=entity_id,
            status=ClassificationStatus.ASSIGNED,
            assignments=tuple(assignments),
            total_n_ORFs=len(orf_ids),
            based_on_n_ORFs=based_on_n_orfs,
            orf_results=tuple(orf_results),
        )



