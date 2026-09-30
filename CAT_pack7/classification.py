#!/usr/bin/env python3
"""
Shared taxonomic classification engine for CAT and BAT.
It makes the logic that is currently duplicated between
CAT and BAT centrally accessible.

I copied it over from earlier work in the cat_pack folder

"""
import gzip
import logging
from dataclasses import dataclass
from decimal import Decimal
from enum import Enum, auto
from pathlib import Path
from typing import Mapping, Sequence

from . import tax
from .utils.errors import InputError
from .utils.logging import Status

log = logging.getLogger("CAT_pack")


class ORFStatus(Enum):
    """Classification status of a single predicted ORF.

    ASSIGNED: has taxid and LCA assigned to ORF
    NO_HIT: ORF has no accepted homology hit
    NO_TAXID: ORF accepted homology hits but no usable taxid association
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
    taxid: str | None = None  # TODO: replace with int down the line
    message: str | None = None
    top_bitscore: Decimal | None = None
    lineage: tuple[str, ...] = ()  # TODO: update on lineage internal change


@dataclass(frozen=True, slots=True)
class TaxonomicAssignment:
    """One taxonomic assignment supported for an entity.

    If f < 0.5 this can cause multiple classifications for one contig.
    """
    taxid: str
    lineage: tuple[str, ...]
    lineage_scores: tuple[Decimal, ...]


class TaxonomyNamespace(Enum):
    """Support for multiple taxonomy namespaces in the future."""
    NCBI = auto()


@dataclass(slots=True)
class ClassificationResult:
    """A unified classification result for CAT and BAT."""
    entity_id: str
    status: ClassificationStatus
    taxonomy_namespace: TaxonomyNamespace = TaxonomyNamespace.NCBI
    assignments: tuple[TaxonomicAssignment, ...] = ()
    total_n_ORFs: int = 0
    based_on_n_ORFs: int = 0
    orf_results: tuple[ORFClassification, ...] = ()


class ClassificationEngine:
    """
    copied over the ClassificationEngine i made in the CAT_pack that shares
    the CAT and BAT logic
    """

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

    def classify_orf(
        self, orf_id: str, hits: Sequence[tuple[str, Decimal]] | None
    ) -> ORFClassification:
        if not hits:
            return ORFClassification(orf_id=orf_id, status=ORFStatus.NO_HIT)

        taxid, top_bitscore = tax.find_LCA_for_ORF(
            hits, self.fastaid2taxid, self.taxid2parent)

        if taxid.startswith("no taxid found"):
            return ORFClassification(
                orf_id=orf_id,
                status=ORFStatus.NO_TAXID,
                n_hits=len(hits),
                top_bitscore=top_bitscore,
                message=taxid,
            )

        lineage = tax.find_lineage(taxid, self.taxid2parent)

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
            lca_ORFs, self.taxid2parent, self.fraction)

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


def import_contig_names(fasta_file: Path) -> set[str]:
    log.info(f"Importing contig names from {fasta_file}.")

    contig_names = set()

    with fasta_file.open("r", encoding="utf-8") as f1:
        for line in f1:
            if line.startswith(">"):
                contig = line.split()[0].lstrip(">").rstrip()

                if contig in contig_names:
                    raise InputError(
                        f"your fasta file contains duplicate headers (the "
                        f"part before the first space in the >line). The "
                        f"first duplicate encountered is {contig}, but there "
                        f"might be more...",
                        path=fasta_file,
                    )

                contig_names.add(contig)

    return contig_names


def import_ORFs(proteins_fasta: Path) -> dict[str, list[str]]:
    log.info(f"Parsing ORF file {proteins_fasta}.")

    contig2ORFs: dict[str, list[str]] = {}

    with proteins_fasta.open("r") as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith(">"):
                ORF = line.split()[0].lstrip(">")
                contig = ORF.rsplit("_", 1)[0]

                if contig not in contig2ORFs:
                    contig2ORFs[contig] = []

                contig2ORFs[contig].append(ORF)

    return contig2ORFs



# I did rewrite this
def parse_alignment(alignment_file: Path, one_minus_r: Decimal) -> tuple[dict[str, list[tuple[str, Decimal]]], set[str]]:
    log.info(f"Parsing alignment file {alignment_file}.")

    opener = gzip.open if alignment_file.suffix == ".gz" else open
    ORF2hits: dict[str, list[tuple[str, Decimal]]] = {}
    all_hits: set[str] = set()

    ORF = "first ORF"
    ORF_done = False
    with opener(alignment_file, "rt", encoding='utf-8') as f1:
        for line in f1:
            line = line.rstrip()
            if not line:
                continue

            fields = line.split("\t")


            if fields[0] == ORF and ORF_done:
                # The ORF has already surpassed its minimum allowed bit-score.
                continue

            if fields[0] != ORF:
                # A new ORF is reached.
                ORF = fields[0]
                top_bitscore = Decimal(fields[11])
                ORF2hits[ORF] = []
                ORF_done = False

            bitscore = Decimal(fields[11])
            if bitscore >= one_minus_r * top_bitscore:
                # The hit has a high enough bit-score to be included.
                hit = fields[1]

                ORF2hits[ORF].append((hit, bitscore))
                all_hits.add(hit)
            else:
                # The hit is not included because its bit-score is too low.
                ORF_done = True

    return ORF2hits, all_hits


def format_lineage(lineage: Sequence[str], branches: set[str]) -> str:
    lineage = list(lineage)
    lineage = tax.star_lineage(lineage, branches)
    return ";".join(lineage[::-1])


def check_orfs_match_contigs(contig_names: set[str], contig2ORFs: Mapping[str, Sequence[str]], path: Path) -> None:
    overlap = len(contig_names & set(contig2ORFs))
    if overlap == 0:
        example = next(iter(contig2ORFs.values()))[0] if contig2ORFs else "contig_name_1"
        raise InputError(
            f"no ORFs found that can be traced back to one of the contigs "
            f"in the contigs fasta file: {example}. ORFs should be named "
            f"contig_name_#.",
            path=path,
        )

    rel_overlap = overlap / len(contig_names)
    log.info(
        f"ORFs found on {overlap:,d} / {len(contig_names):,d} contigs "
        f"({rel_overlap * 100:.2f}%)."
    )
    if rel_overlap < 0.97:
        log.warning(
            f"only {rel_overlap * 100:.2f}% contigs found with ORF predictions. This may "
            f"indicate that some contigs were missing from the protein "
            f"prediction. Please make sure that the protein prediction was "
            f"based on all contigs."
        )


def contig_classification(settings, files, report) -> None:
    contig_names = import_contig_names(files.contigs)
    contig2ORFs = import_ORFs(files.proteins_fasta)
    check_orfs_match_contigs(contig_names, contig2ORFs, files.proteins_fasta)

    one_minus_r = (Decimal("100") - settings.range_) / Decimal("100")
    ORF2hits, all_hits = parse_alignment(files.alignment, one_minus_r)

    taxid2parent, taxid2rank = tax.import_nodes(files.nodes)
    fastaid2LCAtaxid = tax.import_fastaid2LCAtaxid(files.fastaid2LCA, all_hits)
    branches = tax.import_taxids_with_multiple_offspring(files.branches)

    log.info(f"CAT is spinning! Files {files.contig_report} and {files.orf_report} are created.")

    engine = ClassificationEngine(
        taxid2parent=taxid2parent,
        fastaid2taxid=fastaid2LCAtaxid,
        fraction=settings.fraction,
    )

    n_classified_contigs = 0
    n_contigs = len(contig_names)
    report("Classify", Status.RUNNING, 0, n_contigs)

    with (
        files.contig_report.open("w") as contig_out,
        files.orf_report.open("w") as orf_out,
    ):
        contig_out.write(
            f"# contig\tclassification\treason\tlineage\t"
            f"lineage scores (f: {float(settings.fraction)})\n"
        )
        orf_out.write(
            f"# ORF\tnumber of hits (r: {settings.range_})\tlineage\ttop bit-score\n"
        )

        for index, contig in enumerate(sorted(contig_names), start=1):
            result = engine.classify_group(
                entity_id=contig,
                orf_ids=contig2ORFs.get(contig, ()),
                orf2hits=ORF2hits,
            )

            for orf_result in result.orf_results:
                if orf_result.status == ORFStatus.NO_HIT:
                    orf_out.write(
                        f"{orf_result.orf_id}\tORF has no hit to database\n"
                    )
                    continue

                if orf_result.status == ORFStatus.NO_TAXID:
                    orf_out.write(
                        f"{orf_result.orf_id}\t{orf_result.n_hits}\t"
                        f"{orf_result.message}\t{orf_result.top_bitscore}\n"
                    )
                    continue

                orf_out.write(
                    f"{orf_result.orf_id}\t{orf_result.n_hits}\t"
                    f"{format_lineage(orf_result.lineage, branches)}\t"
                    f"{orf_result.top_bitscore}\n"
                )


            match result.status:
                case ClassificationStatus.NO_ORFS:
                    contig_out.write(f"{contig}\tno taxid assigned\tno ORFs found\n")
                case ClassificationStatus.NO_HITS:
                    contig_out.write(f"{contig}\tno taxid assigned\tno hits to database\n")
                case ClassificationStatus.NO_TAXIDS:
                    contig_out.write(f"{contig}\tno taxid assigned\thits not found in taxonomy files\n")
                case ClassificationStatus.NO_LINEAGE_SUPPORT:
                    contig_out.write(f"{contig}\tno taxid assigned\tno lineage reached minimum bit-score support\n")
                case _:
                    n_classified_contigs += 1
                    n_assignments = len(result.assignments)
                    for i, assignment in enumerate(result.assignments):
                        scores = [f"{score:.2f}" for score in assignment.lineage_scores]
                        if n_assignments == 1:
                            label = "taxid assigned"
                        else:
                            f"taxid assigned ({i + 1}/{n_assignments})"
                        contig_out.write(
                            f"{contig}\t{label}\t"
                            f"based on {result.based_on_n_ORFs}/{result.total_n_ORFs} ORFs\t"
                            f"{format_lineage(assignment.lineage, branches)}\t"
                            f"{';'.join(scores[::-1])}\n"
                        )

            report("Classify", Status.RUNNING, index, n_contigs)

    percent = n_classified_contigs / n_contigs * 100 if n_contigs else 0
    log.info(
        f"CAT is done! {n_classified_contigs:,d}/{n_contigs:,d} contigs "
        f"({percent:.2f}%) have taxonomy assigned."
    )
    if settings.fraction < Decimal("0.5"):
        log.warning("since f is set to smaller than 0.5, one contig may have "
                    "multiple classifications.")

    # TODO: show results in the table (or in a panel together with the file table)
    #return n_classified_contigs, n_contigs, percent
