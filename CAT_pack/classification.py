"""Shared taxonomic classification engine for CAT and BAT"""
import logging
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from itertools import chain
from typing import Mapping, Sequence

import numpy as np

from . import tax
from .io.parsers import AlignmentInput
from .results import Classification, ClassificationStatus

log = logging.getLogger("CAT_pack")


# DEVNOTE: Add a flowchart in devdocs to make this more clear for non
# Numpy users
class ClassificationEngine:
    """Classify ORFs and combine their total support for one contig or bin."""

    def __init__(
        self, taxonomy: tax.Taxonomy, fastaid2taxid: np.ndarray,
        fraction: Decimal, threads: int = 1,
    ) -> None:
        self.taxonomy = taxonomy
        self.hit2taxon = fastaid2taxid
        self.fraction = fraction
        self.threads = threads

    def classify(
        self, entity_type: str, entity2ORFs: Mapping[str, Sequence[str]], alignment: AlignmentInput,
    ) -> Classification:
        log.info(f"Starting classification calculation: {len(entity2ORFs)} {entity_type}.")
        orf_entities_names = sorted(entity2ORFs)
        n_entities = len(orf_entities_names)

        # Flatten ORFs, ordered by contig/bin.
        ORF_names = list(chain.from_iterable(entity2ORFs[name] for name in orf_entities_names))
        n_ORFs = np.array([len(entity2ORFs[name]) for name in orf_entities_names], dtype=np.int64)
        ORF_starts = np.concatenate(([0], np.cumsum(n_ORFs)))
        ORF2contig_or_bin = np.repeat(np.arange(n_entities), n_ORFs)
        ORF2alignment_index, ORF_has_hits, LCA_per_ORF, n_hits_per_ORF, top_bitscore_per_ORF = (self.per_ORF_stats(ORF_names, alignment))

        # Combine the ORF lowest common ancestors per contig (or bin)
        # and weighted by their tob bit score(s)
        weighted_LCAs = tax.WeightedLCA.find(
            self.taxonomy, ORF2contig_or_bin, LCA_per_ORF,
            top_bitscore_per_ORF, self.fraction, n_entities)
        log.info("Finished Weighted LCA")
        assignment_starts = np.searchsorted(weighted_LCAs.assignment_entity, np.arange(n_entities + 1))

        def count_per_entity(ORF_mask: np.ndarray) -> np.ndarray:
            return np.bincount(ORF2contig_or_bin[ORF_mask], minlength=n_entities)

        n_ORFs_with_hits = count_per_entity(ORF_has_hits)
        n_ORFs_with_taxid = count_per_entity(LCA_per_ORF != tax.ORF.NO_TAXID)
        n_assignments = np.diff(assignment_starts)

        # DEVNOTE: For in the devdocs. np.select picks the first matching condition
        classification_status = np.select(
            [n_ORFs == 0, n_ORFs_with_hits == 0, n_ORFs_with_taxid == 0, n_assignments == 0],
            [status.value for status in (
                ClassificationStatus.NO_ORFS, ClassificationStatus.NO_HITS,
                ClassificationStatus.NO_TAXIDS, ClassificationStatus.NO_LINEAGE_SUPPORT)],
            default=ClassificationStatus.ASSIGNED.value,
        )

        return Classification(
            entity_type=entity_type,
            entities=orf_entities_names,
            status=classification_status,
            orf_offsets=ORF_starts,
            orf_names=ORF_names,
            orf_row=ORF2alignment_index,
            orf_lca=LCA_per_ORF,
            orf_n_hits=n_hits_per_ORF,
            orf_top_bitscore=top_bitscore_per_ORF,
            n_orfs_with_taxid=n_ORFs_with_taxid,
            assignment_offsets=assignment_starts,
            weighted=weighted_LCAs,
            alignment=alignment,
            taxonomy=self.taxonomy,
        )

    def per_ORF_stats(self, ORF_names: Sequence[str], alignment: AlignmentInput):
        n_ORFs = len(ORF_names)

        ORF2alignment_index = np.fromiter(
            (alignment.orf_index.get(ORF, tax.ORF.NO_HIT) for ORF in ORF_names), dtype=np.int64, count=n_ORFs)
        ORF_has_hits = ORF2alignment_index != tax.ORF.NO_HIT

        alignment_indices = ORF2alignment_index[ORF_has_hits]
        hits_start = alignment.orf_offsets[alignment_indices]
        hits_end = alignment.orf_offsets[alignment_indices + 1]

        hit_taxids = self.hit2taxon[alignment.hit_ids]
        log.info(f"hit taxonomy lookup: {len(hit_taxids)}accepted hits")
        n_alignment_ORFs = len(alignment.orf_names)
        workers = min(self.threads, max(1, n_alignment_ORFs // 20_000))

        log.info(f"Starting ORF LCA calculation: with {n_alignment_ORFs} alignment ORFs, "
                f"{len(hit_taxids)} accepted hits, and {workers} workers.")
        if workers > 1:
            ORF_bounds = np.linspace(
                0, n_alignment_ORFs, workers + 1, dtype=np.int64
            )
            hit_bounds = alignment.orf_offsets[ORF_bounds]
            offset_batches = (
                alignment.orf_offsets[start:end + 1] - hit_start
                for start, end, hit_start in zip(
                    ORF_bounds[:-1], ORF_bounds[1:], hit_bounds[:-1]
                )
            )
            hit_batches = (
                hit_taxids[start:end]
                for start, end in zip(hit_bounds[:-1], hit_bounds[1:])
            )
            with ThreadPoolExecutor(max_workers=workers) as pool:
                alignment_ORF_LCAs = np.concatenate(list(pool.map(
                    self.taxonomy.find_LCA_for_ORF,
                    offset_batches, hit_batches)))
        else:
            alignment_ORF_LCAs = self.taxonomy.find_LCA_for_ORF(
                alignment.orf_offsets, hit_taxids
            )

        LCA_per_ORF = np.full(n_ORFs, tax.ORF.NO_TAXID, dtype=np.int32)
        LCA_per_ORF[ORF_has_hits] = alignment_ORF_LCAs[alignment_indices]

        n_hits_per_ORF = np.zeros(n_ORFs, dtype=np.int64)
        n_hits_per_ORF[ORF_has_hits] = hits_end - hits_start

        top_bitscore_per_ORF = np.zeros(n_ORFs, dtype=np.int64)
        top_bitscore_per_ORF[ORF_has_hits] = alignment.bitscores[hits_start]

        return ORF2alignment_index, ORF_has_hits, LCA_per_ORF, n_hits_per_ORF, top_bitscore_per_ORF

