"""Write the shared CAT/BAT classification report format."""
from decimal import Decimal
from typing import TextIO

from CAT_pack.io.parsers import format_bitscore
from CAT_pack.results import Classification, ClassificationStatus, Lineage
from CAT_pack.tax import Taxonomy, ORF
from CAT_pack.utils.typing import EntityType


class ClassificationWriter:
    """Write headers on creation, then write the array results."""

    def __init__(
        self, classification_out: TextIO, orf_out: TextIO | None = None, *,
        entity_type: EntityType, fraction: Decimal, range_: Decimal,
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
        self._lineage_taxonomy: Taxonomy | None = None
        self._lineages: dict[int, str] = {}
        self._write_headers()

    def _write_headers(self) -> None:
        self.classification_out.write(
            f"# {self.entity_type}\tclassification\treason\tlineage\t"
            f"lineage scores (f: {float(self.fraction)})\n"
        )
        if self.orf_out:
            bin_column = "bin\t" if self.is_bin else ""
            self.orf_out.write(
                f"# ORF\t{bin_column}number of hits (r: {self.range_})\tlineage\ttop bit-score\n"
            )

    def write(self, result: Classification, entity: int) -> None:
        entity_id = result.entities[entity]
        orfs = result.get_orfs(entity)
        if self.orf_out is not None:
            for orf in range(orfs.start, orfs.stop):
                self.write_orf(result, orf, entity)

        status = result.get_status(entity)
        if status != ClassificationStatus.ASSIGNED:
            self.classification_out.write(
                f"{entity_id}\tno taxid assigned\t{status.reason}\n"
            )
            return
        if self.entity_type == EntityType.SEQUENCE:
            reason = f"based on {result.n_w_tax[orfs.start]} hits"
        else:
            reason = (
                f"based on {result.n_orfs_with_taxid[entity]}/"
                f"{orfs.stop - orfs.start} ORFs"
            )

        assignments = result.get_assignments(entity)
        weighted = result.weighted
        total = Decimal(int(weighted.total[entity]))
        score_strings: dict[int, str] = {}

        for index, assignment in enumerate(assignments, start=1):
            label = "taxid assigned"
            if len(assignments) > 1:
                label += f" ({index}/{len(assignments)})"

            start, end = weighted.lineage_offsets[assignment:assignment + 2]
            lineage = self._lineage_for(result.taxonomy, int(weighted.assignment_taxon[assignment]))
            score_parts = []
            for value in weighted.lineage_bitscores[start:end][::-1]:
                score = int(value)
                text = score_strings.get(score)
                if text is None:
                    text = f"{Decimal(score) / total:.2f}"
                    score_strings[score] = text
                score_parts.append(text)
            scores = ";".join(score_parts)
            self.classification_out.write(
                f"{entity_id}\t{label}\t{reason}\t{lineage}\t{scores}\n"
            )

    def write_orf(self, result: Classification, orf: int, entity: int) -> None:
        prefix = result.orf_names[orf]
        if self.is_bin:
            prefix += f"\t{result.entities[entity]}"

        row = int(result.orf_row[orf])
        if row == ORF.NO_HIT:
            self.orf_out.write(f"{prefix}\tORF has no hit to database\n")
            return

        taxid = int(result.orf_lca[orf])
        if taxid == ORF.NO_TAXID:
            hits = ";".join(result.alignment.hit_names_of(row))
            lineage = f"no taxid found ({hits})"
        else:
            lineage = self._lineage_for(result.taxonomy, taxid)

        bitscore = format_bitscore(int(result.orf_top_bitscore[orf]))
        self.orf_out.write(
            f"{prefix}\t{result.orf_n_hits[orf]}\t{lineage}\t{bitscore}\n"
        )

    def format_lineage(self, lineage: tuple[str, ...]) -> str:
        return Lineage(lineage).format(self.branches, no_stars=self.no_stars)

    def _lineage_for(self, taxonomy: Taxonomy, taxid: int) -> str:
        if taxonomy is not self._lineage_taxonomy:
            self._lineages.clear()
            self._lineage_taxonomy = taxonomy
        lineage = self._lineages.get(taxid)
        if lineage is None:
            taxids = tuple(taxonomy.taxids[position] for position in taxonomy.find_lineage(taxid))
            lineage = self.format_lineage(taxids)
            self._lineages[taxid] = lineage
        return lineage
