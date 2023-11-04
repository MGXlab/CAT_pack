#!/usr/bin/env python3

import argparse
import datetime
import decimal
import gzip
import multiprocessing
import os
import pathlib
import subprocess
import sys

import check


class PathAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        path = os.path.expanduser(values.rstrip("/"))

        if not path.startswith("/") and not path.startswith("."):
            path = "./{0}".format(path)

        if os.path.isdir(path):
            path = "{0}/".format(path)

        setattr(namespace, self.dest, path)


class DecimalAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, decimal.Decimal(values))


class SuffixAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        bin_suffix = ".{0}".format(values.lstrip("."))

        setattr(namespace, self.dest, bin_suffix)


def timestamp():
    now = datetime.datetime.now()
    str_ = "[{0}]".format(now.strftime("%Y-%m-%d %H:%M:%S"))

    return str_


def add_argument(argument_group, dest, required, default=None, help_=None):
    if dest == "contigs_fasta":
        if help_ is None:
            help_ = "Path to contigs fasta file."
        argument_group.add_argument(
            "-c",
            "--contigs_fasta",
            dest="contigs_fasta",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "db_fasta":
        if help_ is None:
            help_ = "Path to fasta file containing all sequences."
        argument_group.add_argument(
            "--db_fasta",
            dest="db_fasta",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "bin_fasta_or_folder":
        if help_ is None:
            help_ = "Path to bin fasta file or to directory containing bins."
        argument_group.add_argument(
            "-b",
            "--bin_fasta",
            "--bin_folder",
            dest="bin_fasta_or_folder",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "db_dir":
        if help_ is None:
            help_ = ("Path to directory where CAT/BAT database files will "
                    "be created.")
        argument_group.add_argument(
            "--db_dir",
            dest="db_dir",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "database_folder":
        if help_ is None:
            help_ = "Path to directory that contains database files."
        argument_group.add_argument(
            "-d",
            "--database_folder",
            dest="database_folder",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "taxonomy_folder":
        if help_ is None:
            help_ = "Path to directory that contains taxonomy files."
        argument_group.add_argument(
            "-t",
            "--taxonomy_folder",
            dest="taxonomy_folder",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "names_dmp":
        if help_ is None:
            help_ = "Path to names.dmp"
        argument_group.add_argument(
            "--names",
            dest="names_dmp",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "nodes_dmp":
        if help_ is None:
            help_ = "Path to nodes.dmp"
        argument_group.add_argument(
            "--nodes",
            dest="nodes_dmp",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "acc2tax":
        if help_ is None:
            help_ = "Path to accession2taxid.txt file. Can be gzipped."
        argument_group.add_argument(
            "--acc2tax",
            dest="acc2tax",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "cleanup":
        if help_ is None:
            help_ = ("Remove unnecessary files after all data have been "
                     "processed.")
        argument_group.add_argument(
            "--cleanup",
            dest="cleanup",
            required=required,
            action="store_true",
            help=help_
        )
    elif dest == "bin_suffix":
        if help_ is None:
            help_ = "Suffix of bins in bin directory (default: {0}).".format(
                default)
        argument_group.add_argument(
            "-s",
            "--bin_suffix",
            dest="bin_suffix",
            metavar="",
            required=required,
            type=str,
            default=default,
            help=help_,
        )
    elif dest == "r":
        if help_ is None:
            help_ = "r parameter [0-100] (default: {0:.0f}).".format(default)
        argument_group.add_argument(
            "-r",
            "--range",
            dest="r",
            metavar="",
            required=required,
            type=float,
            choices=[i for i in range(101)],
            action=DecimalAction,
            default=default,
            help=help_,
        )
    elif dest == "f":
        if help_ is None:
            help_ = "f parameter [0-0.99] (default: {0:.2f})." "".format(
                default)
        argument_group.add_argument(
            "-f",
            "--fraction",
            dest="f",
            metavar="",
            required=required,
            type=float,
            choices=[i / 100 for i in range(0, 100)],
            action=DecimalAction,
            default=default,
            help=help_,
        )
    elif dest == "out_prefix":
        if help_ is None:
            help_ = "Prefix for output files (default: {0}).".format(default)
        argument_group.add_argument(
            "-o",
            "--out_prefix",
            dest="out_prefix",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "db":
        if help_ is None:
            help_ = "Either nr  or GTDB."
        argument_group.add_argument(
            "--db",
            dest="db",
            metavar="",
            required=required,
            type=str,
            choices=["nr", "GTDB"],
            default=None,
            help=help_
        )
    elif dest == "output_dir":
        if help_ is None:
            help_ = "Path to directory where data will be stored."
        argument_group.add_argument(
            "-o",
            "--output_dir",
            dest="output_dir",
            metavar="",
            required=required,
            type=lambda p: pathlib.Path(p).resolve(),
            help=help_
        )
    elif dest == "proteins_fasta":
        if help_ is None:
            help_ = (
                "Path to predicted proteins fasta file. If supplied, the "
                "protein prediction step is skipped."
            )
        argument_group.add_argument(
            "-p",
            "--proteins_fasta",
            dest="proteins_fasta",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "alignment_file":
        if help_ is None:
            help_ = (
                "Path to alignment table. If supplied, the alignment step is "
                "skipped and classification is carried out directly. A "
                "predicted proteins fasta file should also be supplied with "
                "argument [-p / --proteins]."
            )
        argument_group.add_argument(
            "-a",
            "--diamond_alignment",
            dest="alignment_file",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "common_prefix":
        if help_ is None:
            help_ = "Prefix for all files that will be created"
        argument_group.add_argument(
            "--common_prefix",
            dest="common_prefix",
            metavar="",
            required=required,
            type=str,
            default=default,
            help=help_,
        )
    elif dest == "path_to_prodigal":
        if help_ is None:
            help_ = (
                "Path to Prodigal binaries. Supply if CAT/BAT cannot find "
                "Prodigal"
            )
        argument_group.add_argument(
            "--path_to_prodigal",
            dest="path_to_prodigal",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "path_to_diamond":
        if help_ is None:
            help_ = (
                "Path to DIAMOND binaries. Supply if CAT/BAT cannot find "
                "DIAMOND."
            )
        argument_group.add_argument(
            "--path_to_diamond",
            dest="path_to_diamond",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            default=default,
            help=help_,
        )
    elif dest == "no_stars":
        if help_ is None:
            help_ = "Suppress marking of suggestive taxonomic assignments."
        argument_group.add_argument(
            "--no_stars",
            dest="no_stars",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "force":
        if help_ is None:
            help_ = "Force overwrite existing files."
        argument_group.add_argument(
            "--force",
            dest="force",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "quiet":
        if help_ is None:
            help_ = "Suppress verbosity."
        argument_group.add_argument(
            "-q",
            "--quiet",
            dest="quiet",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "verbose":
        if help_ is None:
            help_ = "Increase verbosity."
        argument_group.add_argument(
            "--verbose",
            dest="verbose",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "no_log":
        if help_ is None:
            help_ = "Suppress log file."
        argument_group.add_argument(
            "--no_log",
            dest="no_log",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "help":
        if help_ is None:
            help_ = "Show this help message and exit."
        argument_group.add_argument("-h", "--help", action="help", help=help_)
    elif dest == "IkwId":
        if help_ is None:
            help_ = "Flag for experimental features."
        argument_group.add_argument(
            "--I_know_what_Im_doing",
            dest="IkwId",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "input_file":
        if help_ is None:
            help_ = "Path to input file."
        argument_group.add_argument(
            "-i",
            "--input_file",
            dest="input_file",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "output_file":
        if help_ is None:
            help_ = "Path to output file."
        argument_group.add_argument(
            "-o",
            "--output_file",
            dest="output_file",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "only_official":
        if help_ is None:
            help_ = (
                "Only output official raxonomic ranks (superkingdom, phylum, "
                "class, order, family, genus, species)."
            )
        argument_group.add_argument(
            "--only_official",
            dest="only_official",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "exclude_scores":
        if help_ is None:
            help_ = (
                "Do not include bit-score support scores in the lineage of a "
                "classification output file."
            )
        argument_group.add_argument(
            "--exclude_scores",
            dest="exclude_scores",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "nproc":
        if help_ is None:
            help_ = "Number of cores to deploy by DIAMOND (default: maximum)."
        argument_group.add_argument(
            "-n",
            "--nproc",
            dest="nproc",
            metavar="",
            required=required,
            type=int,
            default=default,
            help=help_,
        )
    elif dest == "sensitive":
        if help_ is None:
            help_ = "Run DIAMOND in sensitive mode (default: not enabled)."
        argument_group.add_argument(
            "--sensitive",
            dest="sensitive",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "no_self_hits":
        if help_ is None:
            help_ = (
                "Do not report identical self hits by DIAMOND (default: "
                "not enabled)."
            )
        argument_group.add_argument(
            "--no_self_hits",
            dest="no_self_hits",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "block_size":
        if help_ is None:
            help_ = (
                "DIAMOND block-size parameter (default: {0}). Lower numbers "
                "will decrease memory and temporary disk space usage."
                "".format(default)
            )
        argument_group.add_argument(
            "--block_size",
            dest="block_size",
            metavar="",
            required=required,
            type=float,
            default=default,
            help=help_,
        )
    elif dest == "index_chunks":
        if help_ is None:
            help_ = (
                "DIAMOND index-chunks parameter (default: {0}). Set to 4 on "
                "low memory machines. The parameter has no effect on "
                "temporary disk space usage.".format(default)
            )
        argument_group.add_argument(
            "--index_chunks",
            dest="index_chunks",
            metavar="",
            required=required,
            type=int,
            default=default,
            help=help_,
        )
    elif dest == "tmpdir":
        if help_ is None:
            help_ = (
                "Directory for temporary DIAMOND files (default: directory to "
                "which output files are written)."
            )
        argument_group.add_argument(
            "--tmpdir",
            dest="tmpdir",
            metavar="",
            required=required,
            type=str,
            action=PathAction,
            help=help_,
        )
    elif dest == "compress":
        if help_ is None:
            help_ = "Compress DIAMOND alignment file (default: not enabled)."
        argument_group.add_argument(
            "--compress",
            dest="compress",
            required=required,
            action="store_true",
            help=help_,
        )
    elif dest == "top":
        if help_ is None:
            help_ = (
                "DIAMOND top parameter [0-100] (default: {0}). Governs hits "
                "within range of best hit that are written to the alignment "
                "file. This is not the [-r / --range] parameter! See "
                "README.md.".format(default)
            )
        argument_group.add_argument(
            "--top",
            dest="top",
            metavar="",
            required=required,
            type=float,
            choices=[i for i in range(101)],
            default=default,
            help=help_,
        )
    else:
        sys.exit("Unknown parser dest {0}.".format(dest))

    return


def add_all_diamond_arguments(argument_group):
    add_argument(
        argument_group, "nproc", False, default=multiprocessing.cpu_count()
    )
    add_argument(argument_group, "sensitive", False)
    add_argument(argument_group, "no_self_hits", False)
    add_argument(argument_group, "block_size", False, default=12.0)
    add_argument(argument_group, "index_chunks", False, default=1)
    add_argument(argument_group, "tmpdir", False)
    add_argument(argument_group, "compress", False)
    add_argument(argument_group, "top", False, default=11)

    return


def expand_arguments(args):
    if "r" in args:
        setattr(args, "one_minus_r", (100 - args.r) / 100)

    if "out_prefix" in args:
        # Check out_prefix as the log file needs to be written to a valid
        # location.
        error = check.check_out_prefix(args.out_prefix, None, args.quiet)
        if error:
            sys.exit(1)

        log_file = "{0}.log".format(args.out_prefix)

        with open(log_file, "w") as outf1:
            pass

        if not args.tmpdir:
            tmpdir = "{0}/".format(args.out_prefix.rsplit("/", 1)[0])

            setattr(args, "tmpdir", tmpdir)
    else:
        log_file = None

    setattr(args, "log_file", log_file)

    if "db_dir" in args:
        database_folder_path = str(
            pathlib.Path(args.db_dir) / pathlib.Path("db"))
        diamond_database_name = "{0}.dmnd".format(args.common_prefix)
        diamond_database_path = str(
            database_folder_path / pathlib.Path(diamond_database_name))

        taxonomy_folder_path = str(
            pathlib.Path(args.db_dir) / pathlib.Path("tax"))
        fastaid2LCAtaxid_fname = "{0}.fastaid2LCAtaxid".format(
            args.common_prefix)
        fastaid2LCAtaxid_path = database_folder_path / pathlib.Path(
            fastaid2LCAtaxid_fname)
        fastaid2LCAtaxid_file = str(fastaid2LCAtaxid_path)

        taxids_with_multiple_offspring_fname = (
            "{0}.taxids_with_multiple_offspring".format(args.common_prefix))
        taxids_with_multiple_offspring_path = (
            database_folder_path
            / pathlib.Path(taxids_with_multiple_offspring_fname)
        )
        taxids_with_multiple_offspring_file = str(
            taxids_with_multiple_offspring_path)

        setattr(args, "database_folder", database_folder_path)
        setattr(args, "taxonomy_folder", taxonomy_folder_path)
        setattr(args, "diamond_database", diamond_database_path)
        setattr(args, "fastaid2LCAtaxid_file", fastaid2LCAtaxid_file)
        setattr(
            args,
            "taxids_with_multiple_offspring_file",
            taxids_with_multiple_offspring_file,
        )

    if "taxonomy_folder" in args and not "db_dir" in args:
        setattr(
            args,
            "taxonomy_folder",
            "{0}/".format(args.taxonomy_folder.rstrip("/")),
        )

        explore_taxonomy_folder(args)
    if "database_folder" in args and not "db_dir" in args:
        setattr(
            args,
            "database_folder",
            "{0}/".format(args.database_folder.rstrip("/")),
        )

        explore_database_folder(args)

    if "bin_fasta_or_folder" in args:
        if os.path.isfile(args.bin_fasta_or_folder):
            setattr(args, "bin_fasta", args.bin_fasta_or_folder)
        else:
            setattr(args, "bin_folder", args.bin_fasta_or_folder)

    return


def explore_taxonomy_folder(args):
    nodes_dmp = None
    names_dmp = None
    prot_accession2taxid_file = None

    if os.path.isdir(args.taxonomy_folder):
        for file_ in os.listdir(args.taxonomy_folder):
            if file_ == "nodes.dmp":
                nodes_dmp = "{0}{1}".format(args.taxonomy_folder, file_)
            elif file_ == "names.dmp":
                names_dmp = "{0}{1}".format(args.taxonomy_folder, file_)
            # No need to check for this.
            elif file_.endswith("prot.accession2taxid.FULL.gz"):
                prot_accession2taxid_file = "{0}{1}".format(
                    args.taxonomy_folder, file_)
            elif (
                file_.endswith("prot.accession2taxid.gz")
                and prot_accession2taxid_file is None
            ):
                # Legacy prot_accession2taxid_file.
                prot_accession2taxid_file = "{0}{1}".format(
                    args.taxonomy_folder, file_)

    setattr(args, "nodes_dmp", nodes_dmp)
    setattr(args, "names_dmp", names_dmp)
    setattr(args, "prot_accession2taxid_file", prot_accession2taxid_file)

    return


def explore_database_folder(args):
    fasta_file = None
    diamond_database = None
    fastaid2LCAtaxid_file = None
    taxids_with_multiple_offspring_file = None

    if os.path.isdir(args.database_folder):
        for file_ in os.listdir(args.database_folder):
            if file_.endswith(
                (".fa", ".fasta", ".fna", "fa.gz", "fasta.gz", "fna.gz")):
                fasta_file = "{0}{1}".format(args.database_folder, file_)
            elif file_.endswith(".dmnd"):
                diamond_database = "{0}{1}".format(args.database_folder, file_)
            elif file_.endswith("fastaid2LCAtaxid"):
                fastaid2LCAtaxid_file = "{0}{1}".format(
                    args.database_folder, file_)
            elif file_.endswith("taxids_with_multiple_offspring"):
                taxids_with_multiple_offspring_file = "{0}{1}".format(
                    args.database_folder, file_)

    setattr(args, "db_fasta", fasta_file)
    setattr(args, "diamond_database", diamond_database)
    setattr(args, "fastaid2LCAtaxid_file", fastaid2LCAtaxid_file)
    setattr(
        args,
        "taxids_with_multiple_offspring_file",
        taxids_with_multiple_offspring_file,
    )

    return


def print_variables(args, step_list=None):
    if args.verbose:
        arguments = [
            "{0}: {1}".format(k, v) for k, v in sorted(vars(args).items())
        ]
        message = (
            "\n-----------------\n\n"
            "Full list of arguments:\n"
            "{0}".format("\n".join(arguments))
        )
        give_user_feedback(message, args.log_file, args.quiet, show_time=False)

        if step_list is not None:
            message = "\nStep list: {0}".format(step_list)
            give_user_feedback(
                message, args.log_file, args.quiet, show_time=False)

        message = "\n-----------------\n"
        give_user_feedback(message, args.log_file, args.quiet, show_time=False)

    return


def give_user_feedback(
    message,
    log_file=None,
    quiet=False,
    show_time=True,
    error=False,
    warning=False
):
    if error:
        message = "ERROR: {0}".format(message)
        
    if warning:
        message = "WARNING: {0}".format(message)

    if show_time:
        message = "{0} {1}".format(timestamp(), message)

    message = "{0}\n".format(message)

    if log_file:
        with open(log_file, "a") as outf1:
            outf1.write(message)

    if not quiet and not error:
        sys.stdout.write(message)

    if not quiet and error:
        sys.stderr.write(message)

    return


def run_prodigal(
    path_to_prodigal,
    contigs_fasta,
    proteins_fasta,
    proteins_gff,
    log_file,
    quiet,
):
    message = (
        "Running Prodigal for ORF prediction. Files {0} and {1} will be "
        "generated. Do not forget to cite Prodigal when using CAT or BAT in "
        "your publication.".format(proteins_fasta, proteins_gff)
    )
    give_user_feedback(message, log_file, quiet)

    try:
        command = [
            path_to_prodigal,
            "-i", contigs_fasta,
            "-a", proteins_fasta,
            "-o", proteins_gff,
            "-p", "meta",
            "-g", "11",
            "-q",
            "-f", "gff"
        ]
        subprocess.check_call(command)
    except:
        message = "Prodigal finished abnormally."
        give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = "ORF prediction done!"
    give_user_feedback(message, log_file, quiet)

    return


def run_diamond(args):
    if args.sensitive:
        mode = "sensitive"
    else:
        mode = "fast"

    if args.compress:
        compression = "1"
    else:
        compression = "0"

    message = (
        "Homology search with DIAMOND is starting. Please be patient. Do not "
        "forget to cite DIAMOND when using CAT or BAT in your publication.\n"
        "\t\t\tquery: {0}\n"
        "\t\t\tdatabase: {1}\n"
        "\t\t\tmode: {2}\n"
        "\t\t\ttop: {3}\n"
        "\t\t\tno-self-hits: {4}\n"
        "\t\t\tnumber of cores: {5}\n"
        "\t\t\tblock-size (billions of letters): {6}\n"
        "\t\t\tindex-chunks: {7}\n"
        "\t\t\ttmpdir: {8}\n"
        "\t\t\tcompress: {9}".format(
            args.proteins_fasta,
            args.diamond_database,
            mode,
            args.top,
            args.no_self_hits,
            args.nproc,
            args.block_size,
            args.index_chunks,
            args.tmpdir,
            compression
        )
    )
    give_user_feedback(message, args.log_file, args.quiet)

    try:
        command = [
            args.path_to_diamond, "blastp",
            "-d", args.diamond_database,
            "-q", args.proteins_fasta,
            "--top", str(args.top),
            "--matrix", "BLOSUM62",
            "--evalue", "0.001",
            "-o", args.alignment_file,
            "-p", str(args.nproc),
            "--block-size", str(args.block_size),
            "--index-chunks", str(args.index_chunks),
            "--tmpdir", args.tmpdir,
            "--compress", compression
        ]

        if not args.verbose:
            command += ["--quiet"]

        if args.sensitive:
            command += ["--sensitive"]

        if args.no_self_hits:
            command += ["--no-self-hits"]

        subprocess.check_call(command)
    except:
        message = "DIAMOND finished abnormally."
        give_user_feedback(message, args.log_file, args.quiet, error=True)

        sys.exit(1)

    if args.compress:
        setattr(args, "alignment_file", "{0}.gz".format(args.alignment_file))

    message = "Homology search done! File {0} created.".format(
        args.alignment_file)
    give_user_feedback(message, args.log_file, args.quiet)

    return


def import_contig_names(fasta_file, log_file, quiet):
    message = "Importing contig names from {0}.".format(fasta_file)
    give_user_feedback(message, log_file, quiet)

    contig_names = set()

    with open(fasta_file, "r") as f1:
        for line in f1:
            if line.startswith(">"):
                contig = line.split()[0].lstrip(">").rstrip()

                if contig in contig_names:
                    message = (
                        "your fasta file contains duplicate headers. The "
                        "first duplicate encountered is {0}, but there might "
                        "be more...".format(contig)
                    )
                    give_user_feedback(message, log_file, quiet, error=True)

                    sys.exit(1)

                contig_names.add(contig)

    return contig_names


def import_ORFs(proteins_fasta, log_file, quiet):
    message = "Parsing ORF file {0}".format(proteins_fasta)
    give_user_feedback(message, log_file, quiet)

    contig2ORFs = {}

    with open(proteins_fasta, "r") as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith(">"):
                ORF = line.split()[0].lstrip(">")
                contig = ORF.rsplit("_", 1)[0]

                if contig not in contig2ORFs:
                    contig2ORFs[contig] = []

                contig2ORFs[contig].append(ORF)

    return contig2ORFs


def parse_tabular_alignment(alignment_file, one_minus_r, log_file, quiet):
    message = "Parsing alignment file {0}.".format(alignment_file)
    give_user_feedback(message, log_file, quiet)

    compressed = False
    if alignment_file.endswith(".gz"):
        compressed = True

        f1 = gzip.open(alignment_file, "rb")
    else:
        f1 = open(alignment_file, "r")

    ORF2hits = {}
    all_hits = set()

    ORF = "first ORF"
    ORF_done = False
    for line in f1:
        if compressed:
            line = line.decode("utf-8")

        if line.startswith(ORF) and ORF_done == True:
            # The ORF has already surpassed its minimum allowed bit-score.
            continue

        line = line.rstrip().split("\t")

        if not line[0] == ORF:
            # A new ORF is reached.
            ORF = line[0]
            top_bitscore = decimal.Decimal(line[11])
            ORF2hits[ORF] = []

            ORF_done = False

        bitscore = decimal.Decimal(line[11])

        if bitscore >= one_minus_r * top_bitscore:
            # The hit has a high enough bit-score to be included.
            hit = line[1]

            ORF2hits[ORF].append(
                (hit, bitscore),)
            all_hits.add(hit)
        else:
            # The hit is not included because its bit-score is too low.
            ORF_done = True

    f1.close()

    return (ORF2hits, all_hits)


def is_gz(file_path):
    """Check if given file_paht is gzipped based on suffix."""
    if isinstance(file_path, pathlib.Path):
        file_path = file_path.name
    return file_path.endswith(".gz") or file_path.endswith(".z")


def optionally_compressed_handle(file_path, mode):
    """Return an appropriate file handle to operate on.

    Arguments:
      file_path: str or PathLike: File path.
      mode: str: The passed mode to open the file on.

    Return:
        A file handle either gzip opened or plainly opened for
        reading/writing/appending in text mode.
    """
    if mode == "r" or mode == "rb":
        mode = "rt"
    if mode == "w" or mode == "wb":
        mode = "wt"
    if mode == "a" or mode == "ab":
        mode = "at"

    if is_gz(file_path):
        return gzip.open(file_path, mode=mode)
    else:
        return open(file_path, mode=mode)


if __name__ == "__main__":
    sys.exit("Run \'CAT\' to run CAT or BAT.")
