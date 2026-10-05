"""Write the shared CAT/BAT classification report format."""
from decimal import Decimal
from typing import TextIO

from .results import ClassificationResult, Lineage, ORFClassification, ORFStatus


class ClassificationWriter:
    """Write headers on creation, then write and format the results"""

    def __init__(
        self, classification_out: TextIO, orf_out: TextIO, *,
        entity_type: str, fraction: Decimal, range_: Decimal,
        branches: set[str], no_stars: bool = False,
    ) -> None:
        self.classification_out = classification_out
        self.orf_out = orf_out
        self.entity_type = entity_type
        self.is_bin = entity_type == "bin"
        self.fraction = fraction
        self.range_ = range_
        self.branches = branches
        self.no_stars = no_stars
        self._write_headers()

    def _write_headers(self) -> None:
        self.classification_out.write(
            f"# {self.entity_type}\tclassification\treason\tlineage\t"
            f"lineage scores (f: {float(self.fraction)})\n"
        )
        bin_column = ""
        if self.is_bin:
            bin_column = "bin\t"
        self.orf_out.write(
            f"# ORF\t{bin_column}number of hits (r: {self.range_})\tlineage\ttop bit-score\n"
        )

    def write(self, result: ClassificationResult) -> None:
        for orf_result in result.orf_results:
            self.write_orf(orf_result, result.entity_id)

        if not result.assigned:
            self.classification_out.write(
                f"{result.entity_id}\tno taxid assigned\t{result.reason}\n"
            )
            return

        n_assignments = len(result.assignments)
        reason = result.reason
        for index, assignment in enumerate(result.assignments, start=1):
            label = "taxid assigned"
            if n_assignments > 1:
                label += f" ({index}/{n_assignments})"
            lineage = self.format_lineage(assignment.lineage)
            formatted_scores = []
            for score in reversed(assignment.lineage_scores):
                formatted_scores.append(f"{score:.2f}")
            scores = ";".join(formatted_scores)
            self.classification_out.write(
                f"{result.entity_id}\t{label}\t{reason}\t{lineage}\t{scores}\n"
            )

    def write_orf(self, result: ORFClassification, entity_id: str) -> None:
        prefix = result.orf_id
        if self.is_bin:
            prefix += f"\t{entity_id}"
        if result.status == ORFStatus.NO_HIT:
            self.orf_out.write(f"{prefix}\tORF has no hit to database\n")
            return
        if result.status == ORFStatus.NO_TAXID:
            lineage = result.message
        else:
            lineage = self.format_lineage(result.lineage)
        self.orf_out.write(f"{prefix}\t{result.n_hits}\t{lineage}\t{result.top_bitscore}\n")

    def format_lineage(self, lineage: tuple[str, ...]) -> str:
        return Lineage(lineage).format(self.branches, no_stars=self.no_stars)
