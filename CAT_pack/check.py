#!/usr/bin/env python3

import hashlib
import os
try:
    import pyrodigal
except ImportError:
    pyrodigal = None
import subprocess
import sys

import shared


def check_md5_gz(gz_file, md5_file, log_file, quiet):
    message = "Checking file integrity via MD5 checksum."
    shared.give_user_feedback(message, log_file, quiet)

    with open(md5_file, "r") as f1:
        md5_exp = f1.read().strip().split(" ")[0]

    if md5_exp == "":
        message = ("no MD5 found in {0}. Integrity of {1} cannot be "
                "established.".format(md5_file, gz_file))
        shared.give_user_feedback(message, log_file, quiet, warning=True)
    else:
        md5 = gz_md5(gz_file)

        if md5 != md5_exp:
            message = "MD5 of {0} does not check out.".format(gz_file)
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
        else:
            message = "MD5 of {0} checks out.".format(gz_file)
            shared.give_user_feedback(message, log_file, quiet)

    return


def gz_md5(input_gz, block_size=4096):
    message = "Calculating md5sum for file {0}".format(input_gz)
    shared.give_user_feedback(message)
    md5 = hashlib.md5()
    with open(input_gz, "rb") as f1:
        for chunk in iter(lambda: f1.read(block_size), b""):
            md5.update(chunk)

    return md5.hexdigest()


def check_memory(Gb):
    total_memory = None
    error = False
    
    if sys.platform == "linux" or sys.platform == "linux2":
        # It's a Linux!
        meminfo_file = "/proc/meminfo"
        with open(meminfo_file, "r") as f1:
            for line in f1:
                if line.startswith("MemTotal:"):
                    mem = int(line.split(" ")[-2])

                    # Mem is given in Kb, convert to Gb.
                    total_memory = mem / 2 ** 20
    elif sys.platform == "darwin":
        # It's a Mac!
        meminfo = subprocess.check_output(["sysctl", "hw.memsize"])
        mem = int(meminfo.decode("utf-8").rstrip().split(" ")[-1])
        
        # Mem is given in b, convert to Gb.
        total_memory = mem / 2 ** 30
        
    if total_memory < Gb:
        error = True
        
    return ("{0:.1f}".format(total_memory), error)


def check_out_prefix(out_prefix, log_file, quiet):
    error = False

    if os.path.isdir(out_prefix):
        message = ("prefix for output files ({0}) is a directory."
                "".format(out_prefix))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    dir_ = out_prefix.rsplit("/", 1)[0]

    if not os.path.isdir(dir_):
        message = ("cannot find output directory {0} to which output files "
                "should be written.".format(dir_))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_pyrodigal_install(log_file, quiet, show_time=True):
    error = False

    try:
        version = pyrodigal.__version__

        message = "Pyrodigal found: pyrodigal v{0}.".format(version)
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time)
    except AttributeError:
        message = ("cannot find Pyrodigal. Please check whether it is "
                "installed.")
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time, error=True)

        error = True

    return error


def check_diamond_binaries(path_to_diamond, log_file, quiet, show_time=True):
    error = False

    try:
        p = subprocess.Popen([path_to_diamond, "--version"],
                stdout=subprocess.PIPE)
        c = p.communicate()
        output = c[0].decode().rstrip()

        message = "DIAMOND found: {0}.".format(output)
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time)
    except OSError:
        message = ("cannot find DIAMOND. Please check whether it is "
                "installed or the path to the binaries is provided.")
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time, error=True)

        error = True

    return error


def check_mmseqs2_binaries(path_do_mmseqs2, log_file, quiet, show_time=True):
    error = False

    try:
        p = subprocess.Popen([path_do_mmseqs2, "-h"], stdout=subprocess.PIPE)
        c = p.communicate()
        output = c[0].decode().split("\n")[5]

        message = "MMseqs2 found: {0}.".format(output)
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time)
    except OSError:
        message = ("cannot find MMseqs2. Please check whether it is "
                "installed.")
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time, error=True)

        error = True

    return error


def check_bwa_binaries(path_to_bwa, log_file, quiet, show_time=True):
    error = False

    try:
        p = subprocess.Popen([path_to_bwa],
                             stderr=subprocess.PIPE)
        for line in p.stderr:
            line=line.decode("utf-8")
            if line.startswith("Version"):
                output = line.rstrip()
                message = "bwa found: {0}.".format(output)
                shared.give_user_feedback(
                        message, log_file, quiet, show_time=show_time)
    except OSError:
        message = ("cannot find bwa. Please check whether it is "
                "installed or the path to the binaries is provided.")
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time, error=True)

        error = True

    return error


def check_samtools_binaries(path_to_samtools, log_file, quiet, show_time=True):
    error = False

    try:
        p = subprocess.Popen([path_to_samtools, "--version"],
                             stdout=subprocess.PIPE)
        c = p.communicate()
        output = c[0].decode().split("\n")[0].rstrip()

        message = "Samtools found: {0}.".format(output)
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time)
    except OSError:
        message = ("cannot find samtools. Please check whether it is "
                "installed or the path to the binaries is provided.")
        shared.give_user_feedback(
                message, log_file, quiet, show_time=show_time, error=True)

        error = True

    return error


def check_bin_folder(bin_folder, bin_suffix, log_file, quiet):
    error = False

    if not os.path.isdir(bin_folder):
        message = "cannot find the bin folder."
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

        return error
    
    tmp = []
    with os.scandir(bin_folder) as it:
        for entry in it:
            if entry.name.startswith("."):
                # Skip hidden files.
                continue

            if not entry.name.endswith(bin_suffix):
                continue

            if ".concatenated." in entry.name:
                # Skip concatenated contig fasta and predicted protein fasta files
                # from earlier runs.
                continue
            
            tmp.append(entry.name)

    if len(tmp) == 0:
        message = (
                "no bins found with suffix {0} in bin folder. You can set the "
                "suffix with the [-s / --bin_suffix] argument.".format(
                    bin_suffix)
                )
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_bin_fasta(bin_fasta, log_file, quiet):
    error = False

    if check_fasta(bin_fasta, log_file, quiet):
        error = True

    return error


def check_folders_for_run(
        taxonomy_folder,
        nodes_dmp,
        names_dmp,
        database_folder,
        aligner,
        diamond_database,
        mmseqs2_database,
        fastaid2LCAtaxid_file,
        taxids_with_multiple_offspring_file,
        step_list,
        log_file,
        quiet
        ):
    error = False

    if not os.path.isdir(taxonomy_folder):
        message = "cannot find the taxonomy folder."
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    else:
        if not nodes_dmp or not names_dmp:
            message = ("nodes.dmp and / or names.dmp not found in the "
                    "taxonomy folder.")
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    if not os.path.isdir(database_folder):
        message = "cannot find the database folder."
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    else:
        if (
                aligner.lower() == "diamond" and
                not diamond_database and
                "align" in step_list
                ):
            message = "DIAMOND database not found in database folder."
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if (
                aligner.lower() == "mmseqs2" and
                not mmseqs2_database and
                "align" in step_list
                ):
            message = "MMseqs2 database not found in database folder."
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if not fastaid2LCAtaxid_file:
            message = "file fastaid2LCAtaxid is not found in database folder."
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if not taxids_with_multiple_offspring_file:
            message = ("file taxids_with_multiple_offspring not found in "
                    "database folder.")
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    return error


def check_output_file(output_file, log_file, quiet):
    error = False

    if os.path.isfile(output_file):
        message = (
                "output file {0} already exists. You can choose to overwrite "
                "existing files with the [--force] argument.".format(
                    output_file)
                )
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_input_file(input_file, log_file, quiet):
    error = False

    if not os.path.isfile(input_file):
        message = "input file {0} does not exist.".format(input_file)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_in_and_output_file(input_file, output_file, log_file, quiet):
    error = False

    if input_file == output_file:
        message = "input file and output file cannot be the same."
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_top(top, r, log_file, quiet):
    error = False
    
    if top <= r:
        message = "[--top] should be higher than [-r / --range]."
        shared.give_user_feedback(message, log_file, quiet, error=True)
        
        error = True

    return error


def check_fasta(file_, log_file, quiet):
    error = False

    if not os.path.isfile(file_):
        error = True
    else:
        with open(file_, "r") as f1:
            for n, line in enumerate(f1):
                if n == 0:
                    if not line.startswith(">"):
                        error = True

                    break

    if error:
        message = "{0} is not a fasta file.".format(file_)
        shared.give_user_feedback(message, log_file, quiet, error=True)

    return error


def check_whether_ORFs_are_based_on_contigs(
        contig_names, contig2ORFs, log_file, quiet):
    overlap = len(contig_names & set(contig2ORFs))

    if overlap == 0:
        message = (
                "no ORFs found that can be traced back to one of the contigs "
                "in the contigs fasta file: {0}. ORFs should be named "
                "contig_name_#.".format(contig2ORFs[contig][0])
                )
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    rel_overlap = overlap / len(contig_names)
    message = "ORFs found on {0:,d} / {1:,d} contigs ({2:.2f}%).".format(
            overlap, len(contig_names), rel_overlap * 100)
    shared.give_user_feedback(message, log_file, quiet)

    if rel_overlap < 0.97:
        message = (
                "only {0:.2f}% contigs found with ORF predictions. This may "
                "indicate that some contigs were missing from the protein "
                "prediction. Please make sure that the protein prediction was "
                "based on all contigs.".format(rel_overlap * 100)
                )
        shared.give_user_feedback(message, log_file, quiet, warning=True)

    return
            
            
if __name__ == "__main__":
    sys.exit("Run \'CAT_pack\' to run CAT, BAT, or RAT.")
