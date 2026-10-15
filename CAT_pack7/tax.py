#!/usr/bin/env python3

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import numpy as np

from .io.reference_mapping import read_fastaid2LCAtaxid
from .utils.errors import InputError

log = logging.getLogger("CAT_pack")

NO_TAXID = -1

class Taxonomy:
    def __init__(self, parents, taxid_index, taxids):
        self.parents = parents
        self.taxid_index = taxid_index
        self.taxids = taxids
        self.no_rank = 1
        self.no_taxid = NO_TAXID
        self.prepare_LCA()

    @classmethod
    def import_nodes(cls, nodes_dmp: Path):
        log.info(f"Loading file {nodes_dmp}.")

        taxid2parent = {}

        with open(nodes_dmp, "r") as f1:
            for line in f1:
                line = line.split("\t")
                taxid2parent[line[0]] = line[2]

        taxids = list(taxid2parent)

        # External taxid to internal position index
        taxid_index = {taxid: index for index, taxid in enumerate(taxids)}

        # internal position to parents internal position
        parents = np.fromiter((taxid_index[taxid2parent[taxid]] for taxid in taxids), dtype=np.int32, count=len(taxids))

        return cls(parents, taxid_index, taxids)

    def prepare_LCA(self):
        parents = self.parents
        taxids_indents = np.arange(len(parents))

        taxid_depths = (parents != taxids_indents).astype(np.int64)
        parents_copy = parents.copy()

        for _ in range(len(parents_copy).bit_length()):
            next_parent = parents_copy[parents_copy]

            if np.array_equal(next_parent, parents_copy):
                break

            taxid_depths = taxid_depths + taxid_depths[parents_copy]
            parents_copy = next_parent

        self.depth = taxid_depths
        self.root = parents_copy


    def taxid_of_hits(self, hits, fastaid2taxid):
        taxids = []

        for hit in hits:
            # Missing mappings use a marker outside the internal taxid positions.
            taxids.append(self.taxid_index.get(fastaid2taxid.get(hit), self.no_taxid))

        return np.array(taxids, dtype=np.int32)

    def find_lineage(self, taxid):
        lineage = []
        taxid = int(taxid)

        while True:
            lineage.append(taxid)

            parent = int(self.parents[taxid])
            if parent == taxid:
                break
            taxid = parent

        return np.array(lineage, dtype=np.int32)

    def find_lca(self, taxids):
        taxids = [int(taxid) for taxid in taxids if taxid != self.no_taxid]

        if not taxids:
            return self.no_taxid

        lineages = [self.find_lineage(taxid) for taxid in taxids]

        overlap = set.intersection(*map(set, lineages))

        for taxid in lineages[0]:
            if taxid in overlap:
                return int(taxid)

        return self.no_taxid


    def find_LCA(self, taxids_a, taxids_b):
        taxids_a = np.asarray(taxids_a, dtype=np.int32)
        taxids_b = np.asarray(taxids_b, dtype=np.int32)


        lca_of_taxids = np.full(taxids_a.shape, self.no_taxid, dtype=np.int32)
        shared_root = self.root[taxids_a] == self.root[taxids_b]

        taxids_a = taxids_a[shared_root]
        taxids_b = taxids_b[shared_root]

        while np.any(taxids_a != taxids_b):
            different_taxids = taxids_a != taxids_b
            depth_difference = (self.depth[taxids_a] - self.depth[taxids_b])
            move_taxids_a = different_taxids & (depth_difference >= 0)
            move_taxids_b = different_taxids & (depth_difference <= 0)

            taxids_a[move_taxids_a] = self.parents[taxids_a[move_taxids_a]]
            taxids_b[move_taxids_b] = self.parents[taxids_b[move_taxids_b]]

        lca_of_taxids[shared_root] = taxids_a
        return lca_of_taxids





    def find_LCA_for_ORF(self, orf_hit_offsets, hits):
        batch_started = perf_counter()
        orf_count = len(orf_hit_offsets) - 1

        if orf_count == 0:
            return np.empty(0, dtype=np.int32)

        orf_of_hit = np.repeat(np.arange(orf_count), np.diff(orf_hit_offsets))

        no_taxid_hit_mask = hits != self.no_taxid
        hits = hits[no_taxid_hit_mask]
        orf_of_hit = orf_of_hit[no_taxid_hit_mask]

        hits_per_orf = np.bincount(orf_of_hit, minlength=orf_count)
        max_hits = int(hits_per_orf.max())
        known_hit_count = len(hits)
        log.info(
            "ORF LCA workload: %d ORFs, %d known hits, max hits/ORF=%d.",
            orf_count, known_hit_count, max_hits,
        )
        hit_start = np.concatenate(([0], np.cumsum(hits_per_orf)))[:-1]
        lca_per_orf = np.full(orf_count, self.no_taxid, dtype=np.int32)

        has_hits = hits_per_orf > 0
        reduction_rounds = 0

        while np.any(hits_per_orf > 1):
            reduced_hits_per_orf = (hits_per_orf + 1) // 2
            reduced_hit_start = np.concatenate(([0], np.cumsum(reduced_hits_per_orf)))[:-1]
            orf_of_reduced_hit = np.repeat(np.arange(orf_count), reduced_hits_per_orf)
            first_hit_pos = hit_start[orf_of_reduced_hit] + 2 * (
                np.arange(len(orf_of_reduced_hit)) - reduced_hit_start[orf_of_reduced_hit]
            )
            reduced_hits = hits[first_hit_pos]
            pair_pos = np.flatnonzero(
                first_hit_pos + 1
                < hit_start[orf_of_reduced_hit] + hits_per_orf[orf_of_reduced_hit]
            )
            taxids_a = reduced_hits[pair_pos]
            taxids_b = hits[first_hit_pos[pair_pos] + 1]
            known_pair = (taxids_a != self.no_taxid) & (taxids_b != self.no_taxid)

            reduced_hits[pair_pos[~known_pair]] = self.no_taxid
            different = known_pair & (taxids_a != taxids_b)
            reduced_hits[pair_pos[different]] = self.find_LCA(
                taxids_a[different], taxids_b[different]
            )

            hits = reduced_hits
            hits_per_orf = reduced_hits_per_orf
            hit_start = reduced_hit_start
            reduction_rounds += 1

        lca_per_orf[has_hits] = hits

        elapsed = perf_counter() - batch_started
        log.info(
            "Classification timing: ORF LCA batch: %.3f s (%d ORFs, %d known hits, max hits/ORF=%d, %d rounds)",
            elapsed, orf_count, known_hit_count, max_hits, reduction_rounds,
            extra={"classification_phase": "orf_lca_batch", "classification_seconds": elapsed},
        )
        return lca_per_orf


@dataclass
class WeightedLCA:
    total: np.ndarray
    assignment_entity: np.ndarray
    assignment_taxon: np.ndarray
    lineage_offsets: np.ndarray
    lineage_taxids: np.ndarray
    lineage_bitscores: np.ndarray

    @classmethod
    def find(
        cls, taxonomy, entity_of_orf, orf_lca,
        bitscore, fraction, n_entities,
    ):
        return find_weighted_LCA(
            taxonomy, entity_of_orf, orf_lca,
            bitscore, fraction, n_entities,
        )

    def lineage_support_of(self, assignment):
        start, end = self.lineage_offsets[assignment:assignment + 2]
        return self.lineage_bitscores[start:end].tolist()


def sum_bitscores(orf_type, lca_per_orf, top_bitscores, n):
    taxid_mask = lca_per_orf != NO_TAXID
    total_bitscore = np.zeros(n, dtype=np.int64)
    np.add.at(total_bitscore, orf_type[taxid_mask], top_bitscores[taxid_mask])
    return total_bitscore

def sum_taxid_bitscores(taxonomy, orf_type, lca_per_orf, top_bitscores):
    log.info("Starting support expansion.")
    stage_started = perf_counter()
    orf_pos = np.flatnonzero(lca_per_orf != NO_TAXID)
    known_orfs = len(orf_pos)

    orf_ids = orf_type[orf_pos]
    taxids = lca_per_orf[orf_pos]
    bitscores = top_bitscores[orf_pos]

    support = [(orf_ids, taxids, bitscores, orf_pos)]

    while True:
        parent_taxids = taxonomy.parents[taxids]
        is_not_root = parent_taxids != taxids

        if not np.any(is_not_root):
            break

        orf_ids = orf_ids[is_not_root]
        taxids = parent_taxids[is_not_root]
        bitscores = bitscores[is_not_root]
        orf_pos = orf_pos[is_not_root]

        support.append((orf_ids, taxids, bitscores, orf_pos))

    orf_ids, taxids, bitscores, orf_pos = (
        np.concatenate(column) for column in zip(*support))
    elapsed = perf_counter() - stage_started
    log.info(
        "Classification timing: support expansion: %.3f s (%d known ORFs, %d expanded rows, %d levels)",
        elapsed, known_orfs, len(taxids), len(support),
        extra={"classification_phase": "support_expansion", "classification_seconds": elapsed},
    )

    log.info("Starting support grouping: %d expanded rows.", len(taxids))
    stage_started = perf_counter()
    pairs, pair_pos = np.unique(
        np.column_stack((orf_ids, taxids)),
        axis=0, return_inverse=True
    )
    elapsed = perf_counter() - stage_started
    log.info(
        "Classification timing: support grouping: %.3f s (%d expanded rows, %d unique pairs)",
        elapsed, len(taxids), len(pairs),
        extra={"classification_phase": "support_grouping", "classification_seconds": elapsed},
    )

    stage_started = perf_counter()
    taxid_bitscores = np.zeros(len(pairs), dtype=np.int64)
    np.add.at(taxid_bitscores, pair_pos, bitscores)

    first_orf = np.full(len(pairs), np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(first_orf, pair_pos, orf_pos)
    elapsed = perf_counter() - stage_started
    log.info(
        "Classification timing: support accumulation: %.3f s (%d unique pairs)",
        elapsed, len(pairs),
        extra={"classification_phase": "support_accumulation", "classification_seconds": elapsed},
    )

    return pairs, taxid_bitscores, first_orf


def find_weighted_LCA(taxonomy, orf_type, lca_per_orf, top_bitscores, f, n):
    stage_started = perf_counter()
    total_bitscore = sum_bitscores(orf_type, lca_per_orf, top_bitscores, n)
    elapsed = perf_counter() - stage_started
    log.info(
        "Classification timing: total bitscores: %.3f s (%d entities)",
        elapsed, n,
        extra={"classification_phase": "total_bitscores", "classification_seconds": elapsed},
    )
    pairs, taxid_bitscores, first_orf = sum_taxid_bitscores(
        taxonomy, orf_type, lca_per_orf, top_bitscores)

    log.info("Starting candidate selection: %d support pairs.", len(pairs))
    stage_started = perf_counter()
    f_numerator, f_denominator = f.as_integer_ratio()

    bitscore_cutoffs = np.fromiter(
        (int(total) * f_numerator // f_denominator
         for total in total_bitscore),
        dtype=np.int64, count=n,
    )

    pair_cutoff = bitscore_cutoffs[pairs[:, 0]]

    canidates = np.flatnonzero(taxid_bitscores > pair_cutoff)

    order = np.lexsort((first_orf[canidates], -taxonomy.depth[pairs[canidates, 1]], pairs[canidates, 0]))

    selected_rows = []
    taxid_trace = set()
    previous_orf_id = None

    for canidate in canidates[order]:
        orf_id, taxid = pairs[canidate]

        if orf_id != previous_orf_id:
            taxid_trace.clear()
            previous_orf_id = orf_id

        if taxid in taxid_trace:
            continue

        selected_rows.append(canidate)
        taxid_trace.update(taxonomy.find_lineage(taxid))

    selected_rows = np.asarray(selected_rows, dtype=np.int64)
    elapsed = perf_counter() - stage_started
    log.info(
        "Classification timing: candidate selection: %.3f s (%d candidates, %d assignments)",
        elapsed, len(canidates), len(selected_rows),
        extra={"classification_phase": "candidate_selection", "classification_seconds": elapsed},
    )

    stage_started = perf_counter()
    assignment_entity = pairs[selected_rows, 0]
    assignment_taxon = pairs[selected_rows, 1]

    lineage_lengths = taxonomy.depth[assignment_taxon] + 1
    lineage_offsets = np.concatenate(([0], np.cumsum(lineage_lengths)))

    lineage_taxids = np.empty(lineage_offsets[-1], dtype=np.int32)
    lineage_bitscores = np.empty(lineage_offsets[-1], dtype=np.int64)

    support_offsets = np.searchsorted(pairs[:, 0], np.arange(n + 1))

    for assignment, (entity, taxid) in enumerate(
        zip(assignment_entity, assignment_taxon)
    ):
        lineage = taxonomy.find_lineage(taxid)
        start, end = lineage_offsets[assignment:assignment + 2]

        support_start, support_end = support_offsets[entity:entity + 2]
        support_rows = support_start + np.searchsorted(
            pairs[support_start:support_end, 1], lineage
        )

        lineage_taxids[start:end] = lineage
        lineage_bitscores[start:end] = taxid_bitscores[support_rows]

    result = WeightedLCA(
        total=total_bitscore,
        assignment_entity=assignment_entity,
        assignment_taxon=assignment_taxon,
        lineage_offsets=lineage_offsets,
        lineage_taxids=lineage_taxids,
        lineage_bitscores=lineage_bitscores,
    )
    elapsed = perf_counter() - stage_started
    log.info(
        "Classification timing: lineage packing: %.3f s (%d assignments, %d lineage entries)",
        elapsed, len(assignment_entity), len(lineage_taxids),
        extra={"classification_phase": "lineage_packing", "classification_seconds": elapsed},
    )
    return result

def convert_to_names(lineage, taxid2rank, taxid2name, scores=None):
    names = []
    for (i, taxid) in enumerate(lineage):
        if "*" in taxid:
            taxid = taxid.rstrip("*")

            starred = True
        else:
            starred = False

        name = taxid2name[taxid]
        rank = taxid2rank[taxid]

        if scores is not None:
            if starred:
                names.append("{0}* ({1}): {2}".format(name, rank, scores[i]))
            else:
                names.append("{0} ({1}): {2}".format(name, rank, scores[i]))
        else:
            if starred:
                names.append("{0}* ({1})".format(name, rank))
            else:
                names.append("{0} ({1})".format(name, rank))
                
    return names


def convert_to_official_names(lineage, taxid2rank, taxid2name, scores=None):
    official_ranks = ["superkingdom", "phylum", "class", "order", "family",
                      "genus", "species"]
    lineage_ranks = [taxid2rank[taxid.rstrip("*")] for taxid in lineage]

    official_names = ["no support"] * 7

    for (i, rank) in enumerate(official_ranks):
        if rank in lineage_ranks:
            index = lineage_ranks.index(rank)

            taxid = lineage[index]

            if "*" in taxid:
                taxid = taxid.rstrip("*")

                starred = True
            else:
                starred = False
                
            name = taxid2name[taxid]

            if scores is not None:
                if starred:
                    official_names[i] = "{0}*: {1}".format(name, scores[index])
                else:
                    official_names[i] = "{0}: {1}".format(name, scores[index])
            else:
                if starred:
                    official_names[i] = "{0}*".format(name)
                else:
                    official_names[i] = name

    # Fill the official lineage with NAs if a lower classification is present.
    index_lowest_classification = 0
    for (i, name) in enumerate(official_names):
        if name != "no support":
            index_lowest_classification = i
            
    for i in range(index_lowest_classification):
        if official_names[i] == "no support":
            official_names[i] = "NA"

    return official_names


def read_nodes(path: Path) -> dict[str, str]:
    log.info(f"Loading file {path}.")
    parents: dict[str, str] = {}
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            fields = line.rstrip("\r\n").split("\t", 3)
            if len(fields) < 3 or fields[1] != "|" or not fields[0] or not fields[2]:
                raise InputError(
                    f"Malformed taxonomy node at line {line_number}.",
                    path=path,
                    hint="Use an NCBI nodes.dmp file with tab-separated taxid and parent fields.",
                )
            taxid, parent = fields[0], fields[2]
            if taxid in parents:
                raise InputError(
                    f"Duplicate taxonomy taxid {taxid} at line {line_number}.",
                    path=path,
                )
            parents[taxid] = parent
    return parents


def import_fastaid2LCAtaxid(fastaid2LCAtaxid_file: Path, all_hits, *, workers=1):
    log.info(f"Loading file {fastaid2LCAtaxid_file}.")
    return read_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits, workers)


def import_taxids_with_multiple_offspring(taxids_with_multiple_offspring_file):
    log.info(f"Loading file {taxids_with_multiple_offspring_file}.")

    with open(taxids_with_multiple_offspring_file, "r", encoding="utf-8") as f1:
        return {taxid for line in f1 if (taxid := line.strip())}


def find_lineage(taxid, parents):
    lineage = []
    taxid = int(taxid)

    while True:
        lineage.append(taxid)

        parent = int(parents[taxid])
        if parent == taxid:
            break
        taxid = parent

    return np.array(lineage, dtype=np.int32)


def find_LCA(list_of_lineages):
    overlap = set.intersection(*map(set, list_of_lineages))

    for taxid in list_of_lineages[0]:
        if taxid in overlap:
            return taxid



def import_names(names_dmp: Path):
    log.info(f"Loading file {names_dmp}.")

    taxid2name = {}

    with open(names_dmp, "r") as f1:
        for line in f1:
            line = line.split("\t")

            if line[6] == "scientific name":
                taxid = line[0]
                name = line[2]

                taxid2name[taxid] = name

    return taxid2name


def find_questionable_taxids(lineage, taxids_with_multiple_offspring):
    questionable_taxids = []

    if lineage == ["1"] or lineage == ["root"]:
        return questionable_taxids

    if len(lineage) == 2 and (lineage[1:] == ["1"] or lineage[1:] == ["root"]):
        return questionable_taxids

    for (i, taxid) in enumerate(lineage):
        taxid_parent = lineage[i + 1]
        if taxid_parent in taxids_with_multiple_offspring:
            return questionable_taxids

        questionable_taxids.append(taxid)


def star_lineage(lineage, taxids_with_multiple_offspring):
    questionable_taxids = find_questionable_taxids(
        lineage, taxids_with_multiple_offspring)

    starred_lineage = [
        taxid if
        taxid not in questionable_taxids else
        "{0}*".format(taxid) for taxid in lineage
    ]

    return starred_lineage

if __name__ == "__main__":
    sys.exit("Run \'CAT_pack\' to run CAT, BAT, or RAT.")
