"""Shared taxonomic classification engine for CAT and BAT"""
from decimal import Decimal
from functools import lru_cache
from typing import Mapping, Sequence

from . import tax
from .results import (
    ClassificationResult, ClassificationStatus, ORFClassification,
    ORFStatus, TaxonomicAssignment
)


class ClassificationEngine:
    """Classify ORFs and combine their total support for one contig or bin."""

    def __init__(
        self,
        *,
        taxid2parent: Mapping[str, str],
        fastaid2taxid: Mapping[str, str],
        fraction: Decimal,
    ) -> None:
        self.taxid2parent = taxid2parent
        self.fastaid2taxid = fastaid2taxid
        self.fraction = fraction
        self.lineage = lru_cache(maxsize=65_536)(
            lambda taxid: tuple(tax.find_lineage(taxid, self.taxid2parent)))

    def classify_orf(
        self, orf_id: str, hits: Sequence[tuple[str, Decimal]] | None
    ) -> ORFClassification:
        if not hits:
            return ORFClassification(orf_id=orf_id, status=ORFStatus.NO_HIT)

        taxid, top_bitscore = tax.find_LCA_for_ORF(
            hits, self.fastaid2taxid, self.taxid2parent, lineage_lookup=self.lineage)

        if taxid.startswith("no taxid found"):
            return ORFClassification(
                orf_id=orf_id,
                status=ORFStatus.NO_TAXID,
                n_hits=len(hits),
                top_bitscore=top_bitscore,
                message=taxid,
            )

        lineage = self.lineage(taxid)

        return ORFClassification(
            orf_id=orf_id,
            status=ORFStatus.ASSIGNED,
            n_hits=len(hits),
            taxid=taxid,
            top_bitscore=top_bitscore,
            lineage=tuple(lineage),
        )

    def classify_group(
        self,
        *,
        entity_id: str,
        orf_ids: Sequence[str],
        orf2hits: Mapping[str, Sequence[tuple[str, Decimal]]],
    ) -> ClassificationResult:
        """Classification of one contig or bin from ORFs."""
        if not orf_ids:
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_ORFS,
            )

        orf_results = []
        lca_ORFs = []

        for orf_id in orf_ids:
            result = self.classify_orf(orf_id, orf2hits.get(orf_id))
            orf_results.append(result)

            if result.status == ORFStatus.NO_HIT:
                continue

            if result.status == ORFStatus.NO_TAXID:
                lca_ORFs.append((result.message, result.top_bitscore))
                continue

            lca_ORFs.append((result.taxid, result.top_bitscore))

        if not lca_ORFs:
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_HITS,
                assignments=(),
                total_n_ORFs=len(orf_ids),
                based_on_n_ORFs=0,
                orf_results=tuple(orf_results),
            )

        lineages, lineages_scores, based_on_n_orfs = tax.find_weighted_LCA(
            lca_ORFs, self.taxid2parent, self.fraction, lineage_lookup=self.lineage)

        if lineages == "no ORFs with taxids found.":
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_TAXIDS,
                assignments=(),
                total_n_ORFs=len(orf_ids),
                based_on_n_ORFs=0,
                orf_results=tuple(orf_results),
            )

        if lineages == "no lineage whitelisted.":
            return ClassificationResult(
                entity_id=entity_id,
                status=ClassificationStatus.NO_LINEAGE_SUPPORT,
                assignments=(),
                total_n_ORFs=len(orf_ids),
                based_on_n_ORFs=0,
                orf_results=tuple(orf_results),
            )

        assignments: list[TaxonomicAssignment] = []
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
