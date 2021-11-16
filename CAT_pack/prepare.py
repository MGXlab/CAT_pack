#!/usr/bin/env python3

import argparse
import datetime
import gzip
import multiprocessing
import os
import pathlib
import shutil
import subprocess
import sys
import tarfile
import urllib.request

import about
import check
import shared
import tax


def parse_arguments():
    date = datetime.datetime.now().strftime("%Y-%m-%d")

    parser = argparse.ArgumentParser(
        prog="CAT prepare",
        description="Construct CAT/BAT database files.",
        usage=(
            "CAT prepare --db_fasta path/to/seqs.fa.gz "
            "--prot2acc path/to/accession2taxid.txt "
            "--names path/to/names.dmp "
            "--nodes path/to/nodes.dmp "
            "[options] [-h / --help]"
        ),
        add_help=False,
    )

    required = parser.add_argument_group("Required arguments")
    shared.add_argument(required, "db_fasta", True)
    shared.add_argument(required, "names_dmp", True)
    shared.add_argument(required, "nodes_dmp", True)
    shared.add_argument(required, "prot2acc", True)

    optional = parser.add_argument_group("Optional arguments")
    shared.add_argument(
        optional,
        "database_folder",
        False,
        default="./CAT_database.{0}".format(date),
        help_=(
            "Name of folder to which database files will be written "
            "(default: CAT_database.{date})"
        ),
    )
    shared.add_argument(optional, "path_to_diamond", False, default="diamond")
    shared.add_argument(optional, 'diamond_database_prefix', False,
                        default="CAT",
                        help_=("Prefix for the .dmnd file to be created "),
                        )
    shared.add_argument(optional, "quiet", False)
    shared.add_argument(optional, "verbose", False)
    shared.add_argument(optional, "no_log", False)
    shared.add_argument(optional, "help", False)

    specific = parser.add_argument_group("DIAMOND specific optional arguments")
    shared.add_argument(
        specific, "nproc", False, default=multiprocessing.cpu_count()
    )

    (args, extra_args) = parser.parse_known_args()

    extra_args = [
        arg for (i, arg) in enumerate(extra_args) if (i, arg) != (0, "prepare")
    ]
    if len(extra_args) > 0:
        sys.exit(
            "error: too much arguments supplied:\n{0}".format(
                "\n".join(extra_args)
            )
        )

    # Add extra arguments.
    setattr(args, "date", date)
    setattr(args, "min_mem", 200)
    shared.expand_arguments(args)

    return args


def memory_bottleneck(args):
    (total_memory, error) = check.check_memory(args.min_mem)
    if error:
        message = (
            "WARNING: "
            "At least {0}GB of memory is recommended for large database "
            "construction (e.g. nr). {1}GB is found on your system. "
            "You can try to find a machine with more memory if you run into "
            "issues or download preconstructed database files from "
            "tbb.bio.uu.nl/bastiaan/CAT_prepare/.".format(
                args.min_mem, total_memory
            )
        )
        shared.give_user_feedback(
            message, args.log_file, args.quiet,
        )

    return


def make_diamond_database(
    path_to_diamond,
    fasta_file,
    diamond_database_prefix,
    nproc,
    log_file,
    quiet,
    verbose,
):
    message = (
        "Constructing DIAMOND database {0}.dmnd from {1} using {2} "
        "cores.".format(diamond_database_prefix, fasta_file, nproc)
    )
    shared.give_user_feedback(message, log_file, quiet)

    command = [
        path_to_diamond,
        "makedb",
        "--in",
        fasta_file,
        "-d",
        diamond_database_prefix,
        "-p",
        str(nproc),
    ]

    if not verbose:
        command += ["--quiet"]

    try:
        subprocess.check_call(command)
    except:
        message = "DIAMOND database could not be created."
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = "DIAMOND database constructed."
    shared.give_user_feedback(message, log_file, quiet)

    return


def import_fasta_headers(fasta, log_file, quiet):
    message = "Loading file {0}.".format(fasta)
    shared.give_user_feedback(message, log_file, quiet)

    fastaid2prot_accessions = {}
    prot_accessions_whitelist = set()

    with shared.optionally_compressed_handle(fasta, "r") as f1:
        for line in f1:

            if not line.startswith(">"):
                continue

            # \x01 == ^A, handles multiple fasta headers
            # Some legacy format for headers
            line = line.lstrip(">").split("\x01")

            prot_accessions = [i.split(" ")[0].strip() for i in line]
            fastaid = prot_accessions[0]

            fastaid2prot_accessions[fastaid] = prot_accessions
            prot_accessions_whitelist.update(prot_accessions)

    return (fastaid2prot_accessions, prot_accessions_whitelist)


def import_prot_accession2taxid(
    prot_accession2taxid_file, prot_accessions_whitelist, log_file, quiet
):
    message = "Loading file {0}.".format(prot_accession2taxid_file)
    shared.give_user_feedback(message, log_file, quiet)

    prot_accession2taxid = {}

    with shared.optionally_compressed_handle(
        prot_accession2taxid_file, "r"
    ) as f1:
        for n, line in enumerate(f1):

            line = line.rstrip().split("\t")

            if n == 0:
                index_1 = line.index("accession.version")
                index_2 = line.index("taxid")

                continue

            prot_accession = line[index_1]

            if prot_accession in prot_accessions_whitelist:
                prot_accession2taxid[prot_accession] = line[index_2]

    return prot_accession2taxid


def make_fastaid2LCAtaxid_file(
    fastaid2LCAtaxid_file,
    fasta_file,
    prot_accession2taxid_file,
    taxid2parent,
    log_file,
    quiet,
):
    (
        fastaid2prot_accessions,
        prot_accessions_whitelist,
    ) = import_fasta_headers(fasta_file, log_file, quiet)
    prot_accession2taxid = import_prot_accession2taxid(
        prot_accession2taxid_file, prot_accessions_whitelist, log_file, quiet
    )

    message = "Finding LCA of all protein accession numbers in fasta headers."
    shared.give_user_feedback(message, log_file, quiet)

    no_taxid = 0
    corrected = 0
    total = 0
    with open(fastaid2LCAtaxid_file, "w") as outf1:
        for fastaid, prot_accessions in fastaid2prot_accessions.items():
            list_of_lineages = []
            for prot_accession in prot_accessions:
                try:
                    taxid = prot_accession2taxid[prot_accession]
                    lineage = tax.find_lineage(taxid, taxid2parent)
                    list_of_lineages.append(lineage)
                except:
                    # This accounts for missing accession numbers in
                    # prot.accession2taxid and missing nodes in nodes.dmp.
                    continue

            total += 1

            if len(list_of_lineages) == 0:
                # This accounts for entries that only contain accession numbers
                # that are missing in prot.accession2taxid or whose taxid is
                # missing in nodes.dmp. NOTE that these entries are thus not
                # present in the output file.
                no_taxid += 1

                continue

            LCAtaxid = tax.find_LCA(list_of_lineages)

            outf1.write("{0}\t{1}\n".format(fastaid, LCAtaxid))

            if (
                fastaid not in prot_accession2taxid
                or LCAtaxid != prot_accession2taxid[fastaid]
            ):
                # If the fastaid cannot be found in prot.accession2taxid, but
                # a taxid is given to the fastaid based on secondary accession
                # numbers, or if the taxid of the header is different from the
                # LCA taxid, it is counted as corrected.
                corrected += 1

    message = (
        "Done! File {0} is created. "
        "{1:,d} of {2:,d} headers ({3:.1f}%) corrected. "
        "{4:,d} headers ({5:.1f}%) do not have a taxid assigned.".format(
            fastaid2LCAtaxid_file,
            corrected,
            total,
            corrected / total * 100,
            no_taxid,
            no_taxid / total * 100,
        )
    )
    shared.give_user_feedback(message, log_file, quiet)

    return


def find_offspring(fastaid2LCAtaxid_file, taxid2parent, log_file, quiet):
    message = "Searching database for taxids with multiple offspring."
    shared.give_user_feedback(message, log_file, quiet)

    taxid2offspring = {}

    with open(fastaid2LCAtaxid_file, "r") as f1:
        for line in f1:
            line = line.rstrip().split("\t")

            taxid = line[1]
            lineage = tax.find_lineage(taxid, taxid2parent)

            for (i, taxid) in enumerate(lineage):
                # The first taxid in the lineage does not have a daughter node.
                if i == 0:
                    continue

                if taxid not in taxid2offspring:
                    taxid2offspring[taxid] = set()

                offspring = lineage[i - 1]

                taxid2offspring[taxid].add(offspring)

    return taxid2offspring


def write_taxids_with_multiple_offspring_file(
    taxids_with_multiple_offspring_file, taxid2offspring, log_file, quiet
):
    message = "Writing {0}.".format(taxids_with_multiple_offspring_file)
    shared.give_user_feedback(message, log_file, quiet)

    with open(taxids_with_multiple_offspring_file, "w") as outf1:
        for taxid in taxid2offspring:
            if len(taxid2offspring[taxid]) >= 2:
                outf1.write("{0}\n".format(taxid))

    return


def prepare(step_list, args):
    shared.print_variables(args, step_list)
    skip_steps = []
    memory_bottleneck(args)

    # This is the root dir
    db_dir = pathlib.Path(args.database_folder).resolve()
    db_dir.mkdir(exist_ok=True)

    # It should contain 
    # 1. a taxonomy folder with names and nodes
    tax_db = db_dir / pathlib.Path("tax")
    tax_db.mkdir(exist_ok=True)

    # Check if names and nodes exist together
    nodes_dmp_fp = tax_db / pathlib.Path("nodes.dmp")
    names_dmp_fp = tax_db / pathlib.Path("names.dmp")
    if not names_dmp_fp.exists() or not nodes_dmp_fp.exists():
        message = ("Copying names.dmp and nodes.dmp in {}".format(tax_db))
        shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=True
        )
        shutil.copyfile(args.names_dmp, names_dmp_fp)
        shutil.copyfile(args.nodes_dmp, nodes_dmp_fp)


    # 2. a dir with the .dmnd and LCA files
    cat_db = db_dir / pathlib.Path("db")

    if cat_db.is_dir():
        if any(cat_db.glob("*.dmnd")):
            message = "A DIAMOND database exists. Skipping creation"
            shared.give_user_feedback(
                message, args.log_file, args.quiet, show_time=False
            )
            step_list.pop(0)
    else:
        message = "Database folder is created at {}".format(cat_db)
        shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=True
        )
        cat_db.mkdir()

    setattr(
        args,
        "diamond_database_prefix",
        cat_db / pathlib.Path(args.diamond_database_prefix)
    )


    if "make_diamond_database" in step_list:
        make_diamond_database(
            args.path_to_diamond,
            args.db_fasta,
            args.diamond_database_prefix,
            args.nproc,
            args.log_file,
            args.quiet,
            args.verbose,
        )

    if (
        "make_fastaid2LCAtaxid_file" in step_list
        or "make_taxids_with_multiple_offspring_file" in step_list
    ):
        taxid2parent, taxid2rank = tax.import_nodes(
            args.nodes_dmp, args.log_file, args.quiet
        )

    if "make_fastaid2LCAtaxid_file" in step_list:
        fname = "{}.fastaid2LCAtaxid".format(args.date)
        fpath = cat_db / pathlib.Path(fname)
        setattr(
            args,
            "fastaid2LCAtaxid_file",
            fpath
        )

        make_fastaid2LCAtaxid_file(
            args.nodes_dmp,
            args.fastaid2LCAtaxid_file,
            args.db_fasta,
            args.prot2acc,
            taxid2parent,
            args.log_file,
            args.quiet,
        )

    if "make_taxids_with_multiple_offspring_file" in step_list:
        fname = "{}.taxids_with_multiple_offspring".format(args.date)
        fpath = cat_db / pathlib.Path(fname)
        setattr(
            args,
            "taxids_with_multiple_offspring_file",
            fpath
            )

        taxid2offspring = find_offspring(
            args.nodes_dmp,
            args.fastaid2LCAtaxid_file,
            taxid2parent,
            args.log_file,
            args.quiet,
        )
        write_taxids_with_multiple_offspring_file(
            args.taxids_with_multiple_offspring_file,
            taxid2offspring,
            args.log_file,
            args.quiet,
        )

    message = "\n-----------------\n\n" "{0} CAT prepare is done!".format(
        shared.timestamp()
    )
    shared.give_user_feedback(
        message, args.log_file, args.quiet, show_time=False
    )

    message = (
        "\nSupply the following arguments to CAT or BAT if you want to "
        "use this database:\n"
        "-d / --database_folder {0}\n"
        "-t / --taxonomy_folder {1}".format(
            cat_db, tax_db
        )
    )
    shared.give_user_feedback(
        message, args.log_file, args.quiet, show_time=False
    )

    return


def run():
    args = parse_arguments()
    step_list = [
        "make_diamond_database",
        "make_fastaid2LCAtaxid_file",
        "make_taxids_with_multiple_offspring_file",
    ]
    prepare(step_list, args)
    return


if __name__ == "__main__":
    sys.exit("Run 'CAT prepare' to construct a CAT/BAT database.")
