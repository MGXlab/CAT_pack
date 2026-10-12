"""Write the shared CAT/BAT classification report format."""
from decimal import Decimal
from typing import TextIO

from .parsers import format_bitscore
from ..results import Classification, ClassificationStatus, Lineage, NO_HIT
from ..tax import NO_TAXID


class ClassificationWriter:
    """Write headers on creation, then write the array results."""

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
        bin_column = "bin\t" if self.is_bin else ""
        self.orf_out.write(
            f"# ORF\t{bin_column}number of hits (r: {self.range_})\tlineage\ttop bit-score\n"
        )

    def write(self, result: Classification, entity: int) -> None:
        entity_id = result.entities[entity]
        orfs = result.get_orfs(entity)
        for orf in range(orfs.start, orfs.stop):
            self.write_orf(result, orf, entity)

        status = result.get_status(entity)
        if status != ClassificationStatus.ASSIGNED:
            self.classification_out.write(
                f"{entity_id}\tno taxid assigned\t{status.reason}\n"
            )
            return

        assignments = result.get_assignments(entity)
        reason = (
            f"based on {result.n_orfs_with_taxid[entity]}/"
            f"{orfs.stop - orfs.start} ORFs"
        )
        weighted = result.weighted
        total = int(weighted.total[entity])

        for index, assignment in enumerate(assignments, start=1):
            label = "taxid assigned"
            if len(assignments) > 1:
                label += f" ({index}/{len(assignments)})"

            start, end = weighted.lineage_offsets[assignment:assignment + 2]
            taxids = tuple(
                result.taxonomy.taxids[taxid]
                for taxid in weighted.lineage_taxids[start:end]
            )
            lineage = self.format_lineage(taxids)
            scores = ";".join(
                f"{Decimal(int(score)) / Decimal(total):.2f}"
                for score in weighted.lineage_bitscores[start:end][::-1]
            )
            self.classification_out.write(
                f"{entity_id}\t{label}\t{reason}\t{lineage}\t{scores}\n"
            )

    def write_orf(self, result: Classification, orf: int, entity: int) -> None:
        prefix = result.orf_names[orf]
        if self.is_bin:
            prefix += f"\t{result.entities[entity]}"

        row = int(result.orf_row[orf])
        if row == NO_HIT:
            self.orf_out.write(f"{prefix}\tORF has no hit to database\n")
            return

        taxid = int(result.orf_lca[orf])
        if taxid == NO_TAXID:
            hits = ";".join(result.alignment.hit_names_of(row))
            lineage = f"no taxid found ({hits})"
        else:
            taxids = tuple(
                result.taxonomy.taxids[position]
                for position in result.taxonomy.find_lineage(taxid)
            )
            lineage = self.format_lineage(taxids)

        bitscore = format_bitscore(int(result.orf_top_bitscore[orf]))
        self.orf_out.write(
            f"{prefix}\t{result.orf_n_hits[orf]}\t{lineage}\t{bitscore}\n"
        )

    def format_lineage(self, lineage: tuple[str, ...]) -> str:
        return Lineage(lineage).format(self.branches, no_stars=self.no_stars)
