#!/usr/bin/env python3

import sys

import shared


def import_nodes(nodes_dmp, log_file, quiet):
    message = "Loading file {0}.".format(nodes_dmp)
    shared.give_user_feedback(message, log_file, quiet)
    
    taxid2parent = {}
    taxid2rank = {}

    with open(nodes_dmp, "r") as f1:
        for line in f1:
            line = line.split("\t")

            taxid = line[0]
            parent = line[2]
            rank = line[4]

            taxid2parent[taxid] = parent
            taxid2rank[taxid] = rank

    return (taxid2parent, taxid2rank)


def import_names(names_dmp, log_file, quiet):
    message = "Loading file {0}.".format(names_dmp)
    shared.give_user_feedback(message, log_file, quiet)

    taxid2name = {}

    with open(names_dmp, "r") as f1:
        for line in f1:
            line = line.split("\t")

            if line[6] == "scientific name":
                taxid = line[0]
                name = line[2]

                taxid2name[taxid] = name

    return taxid2name


def import_fastaid2LCAtaxid(fastaid2LCAtaxid_file, all_hits, log_file, quiet):
    message = "Loading file {0}.".format(fastaid2LCAtaxid_file)
    shared.give_user_feedback(message, log_file, quiet)

    fastaid2LCAtaxid = {}

    with open(fastaid2LCAtaxid_file, "r") as f1:
        for line in f1:
            line = line.rstrip().split("\t")

            if line[0] in all_hits:
                # Only include fastaids that are found in hits.
                fastaid2LCAtaxid[line[0]] = line[1]

    return fastaid2LCAtaxid


def import_taxids_with_multiple_offspring(
    taxids_with_multiple_offspring_file, log_file, quiet):
    message = "Loading file {0}.".format(taxids_with_multiple_offspring_file)
    shared.give_user_feedback(message, log_file, quiet)

    taxids_with_multiple_offspring = set()

    with open(taxids_with_multiple_offspring_file, "r") as f1:
        for line in f1:
            line = line.rstrip()

            taxids_with_multiple_offspring.add(line)

    return taxids_with_multiple_offspring


def find_lineage(taxid, taxid2parent, lineage=None):
    if lineage is None:
        lineage = []

    lineage.append(taxid)

    if taxid2parent[taxid] == taxid:
        return lineage
    else:
        return find_lineage(taxid2parent[taxid], taxid2parent, lineage)
    
    
def find_LCA(list_of_lineages):
    overlap = set.intersection(*map(set, list_of_lineages))

    for taxid in list_of_lineages[0]:
        if taxid in overlap:
            return taxid


def find_LCA_for_ORF(hits, fastaid2LCAtaxid, taxid2parent):
    list_of_lineages = []
    top_bitscore = 0

    for (hit, bitscore) in hits:
        if bitscore > top_bitscore:
            top_bitscore = bitscore
            
        try:
            taxid = fastaid2LCAtaxid[hit]
            lineage = find_lineage(taxid, taxid2parent)

            list_of_lineages.append(lineage)
        except:
            # The fastaid does not have an associated taxid for some reason.
            pass
        
    if len(list_of_lineages) == 0:
        return (
                "no taxid found ({0})".format(";".join([i[0] for i in hits])),
                top_bitscore
                )

    overlap = set.intersection(*map(set, list_of_lineages))

    for taxid in list_of_lineages[0]:
        if taxid in overlap:
            return (taxid, top_bitscore)
        
        
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


def find_weighted_LCA(LCAs_ORFs, taxid2parent, f):
    list_of_lineages = []
    list_of_bitscores = []
    based_on_n_ORFs = 0

    for (taxid, top_bitscore) in LCAs_ORFs:
        if taxid.startswith("no taxid found"):
            # Thus the ORFs that are not classified because they don"t have an
            # associated taxid are not taken into account for the
            # classification of the contig.
            continue
        
        lineage = find_lineage(taxid, taxid2parent)
        
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

    whitelisted_lineages = []
    for taxid in taxid2bitscore:
        if taxid2bitscore[taxid] / sum(list_of_bitscores) > f:
            lineage = find_lineage(taxid, taxid2parent)

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
            longest_lineages.append(whitelisted_lineage)

            scores = [taxid2bitscore[taxid] / sum(list_of_bitscores) for
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


if __name__ == "__main__":
    sys.exit("Run \'CAT_pack\' to run CAT, BAT, or RAT.")
