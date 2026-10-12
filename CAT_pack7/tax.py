#!/usr/bin/env python3

import logging
import sys
from pathlib import Path

import numpy as np

from .io.reference_mapping import read_fastaid2LCAtaxid

log = logging.getLogger("CAT_pack")

class Taxonomy:
    def __init__(self, parents, taxid_index, taxids):
        self.parents = parents
        self.taxid_index = taxid_index
        self.taxids = taxids
        self.no_rank = 1
        self.no_taxid = 0
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
            # Note, no rank is 1 in the nodes.dmp that's why the fallback is 0
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

        lineages = self.find_lineage(taxid for taxid in taxids)

        overlap = set.intersection(*map(set, lineages))

        for taxid in lineages[0]:
            if taxid in overlap:
                return int(taxid)

        return self.no_taxid


    def find_LCA(self, taxids_a, taxids_b):
        taxids_a = np.asarray(taxids_a, dtype=np.int32)
        taxids_b = np.asarray(taxids_b, dtype=np.int32)


        lca_of_taxids = np.full(taxids_a.shape, -1, dtype=np.int32)
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
        orf_count = len(orf_hit_offsets) - 1

        lca = np.full(orf_count, self.no_taxid, dtype=np.int32)

        for orf in range(orf_count):












        
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


def find_weighted_LCA(LCAs_ORFs, taxid2parent, f, *, lineage_lookup=None):
    if lineage_lookup is None:
        lineage_lookup = lambda taxid: find_lineage(taxid, taxid2parent)
    list_of_lineages = []
    list_of_bitscores = []
    based_on_n_ORFs = 0

    for (taxid, top_bitscore) in LCAs_ORFs:
        if taxid.startswith("no taxid found"):
            # Thus the ORFs that are not classified because they don't have an
            # associated taxid are not taken into account for the
            # classification of the contig.
            continue
        
        lineage = lineage_lookup(taxid)
        
        list_of_lineages.append(lineage)
        list_of_bitscores.append(top_bitscore)
        based_on_n_ORFs += 1

    if len(list_of_lineages) == 0:
        return (
                "no ORFs with taxids found.",
                "no ORFs with taxids found.",
                "no ORFs with taxids found."
                )

    taxid2bitscore = {}
    for (i, lineage) in enumerate(list_of_lineages):
        for taxid in lineage:
            if taxid not in taxid2bitscore:
                taxid2bitscore[taxid] = 0

            taxid2bitscore[taxid] += list_of_bitscores[i]

    total_bitscore = sum(list_of_bitscores)
    whitelisted_lineages = []
    for taxid in taxid2bitscore:
        if taxid2bitscore[taxid] / total_bitscore > f:
            lineage = lineage_lookup(taxid)

            whitelisted_lineages.append(lineage)

    if len(whitelisted_lineages) == 0:
        return (
                "no lineage whitelisted.",
                "no lineage whitelisted.",
                "no lineage whitelisted."
                )

    whitelisted_lineages = sorted(whitelisted_lineages,
            key=lambda x: len(x), reverse=True)

    longest_lineages = []
    longest_lineages_scores = []

    taxid_trace = set()
    for whitelisted_lineage in whitelisted_lineages:
        if whitelisted_lineage[0] not in taxid_trace:
            longest_lineages.append(list(whitelisted_lineage))

            scores = [taxid2bitscore[taxid] / total_bitscore for
                      taxid in whitelisted_lineage]
            longest_lineages_scores.append(scores)

            taxid_trace |= set(whitelisted_lineage)

    return (longest_lineages, longest_lineages_scores, based_on_n_ORFs)


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


def import_fastaid2LCAtaxid(fastaid2LCAtaxid_file: Path, all_hits, *, workers=1):
    log.info(f"Loading file {fastaid2LCAtaxid_file}.")
    return read_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits, workers)


def import_taxids_with_multiple_offspring(taxids_with_multiple_offspring_file):
    log.info(f"Loading file {taxids_with_multiple_offspring_file}.")

    taxids_with_multiple_offspring = set()

    with open(taxids_with_multiple_offspring_file, "r") as f1:
        for line in f1:
            line = line.rstrip()

            taxids_with_multiple_offspring.add(line)

    return taxids_with_multiple_offspring


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

if __name__ == "__main__":
    sys.exit("Run \'CAT_pack\' to run CAT, BAT, or RAT.")
