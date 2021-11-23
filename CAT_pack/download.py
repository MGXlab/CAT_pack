import argparse
import datetime
import pathlib
import sys
import tarfile
import urllib.request

import shared
import check

REMOTE_URLS = {
    "nr": [
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz.md5",
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5",
    ],
    "gtdb": [
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122_taxonomy.tsv.gz",
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_taxonomy.tsv.gz",
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM",
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz",
    ],
}


def parse_arguments():
    date = datetime.datetime.now().strftime("%Y-%m-%d")

    parser = argparse.ArgumentParser(
        prog="CAT download",
        description=(
            "Download and preprocess sequence and taxonomy information. "
            "Currently supports data from nr and GTDB"
        ),
        usage="CAT download -db [nr|gtdb] -o output_dir",
        add_help=False,
    )

    required = parser.add_argument_group("Required Arguments")
    shared.add_argument(required, "db", True)
    shared.add_argument(required, "output_dir", True)

    optional = parser.add_argument_group("Otpional Arguments")
    shared.add_argument(optional, "quiet", False)
    shared.add_argument(optional, "verbose", False)
    shared.add_argument(optional, "no_log", False)
    shared.add_argument(optional, "help", False)

    (args, extra_args) = parser.parse_known_args()

    extra_args = [
        arg
        for (i, arg) in enumerate(extra_args)
        if (i, arg) != (0, "download")
    ]
    if len(extra_args) > 0:
        sys.exit(
            "error: too many arguments supplied:\n{0}".format(
                "\n".join(extra_args)
            )
        )

    setattr(args, "date", date)
    shared.expand_arguments(args)
    return args


def download_singleton(target_url, local_path, log_file, quiet):
    """Download a single file to the specified location"""
    try:
        urllib.request.urlretrieve(target_url, local_path)
    except:
        message = "Failed downloading file: {}".format(target_url)
        shared.give_user_feedback(message, log_file, quiet)
        raise


def nr_download(output_dir, date_prefix, log_file, quiet):
    """Download all required nr files in the specified output dir

    All files downloaded are prefixed with a YYYY-MM-DD timestamp
    """

    nr_urls = REMOTE_URLS.get("nr")

    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        existing_files = []
    else:
        existing_files = list([p.resolve() for p in output_dir.iterdir()])

    for url in nr_urls:
        url_leaf = url.split("/")[-1]
        output_basename = "{}.{}".format(date_prefix, url_leaf)
        output_path = output_dir / pathlib.Path(output_basename)
        if output_path in existing_files:
            message = (
                "Skipping download of file {} . It already exists.".format(
                    output_path.name
                )
            )
            shared.give_user_feedback(message, log_file, quiet)
        else:
            message = "Downloading {}".format(url_leaf)
            shared.give_user_feedback(message, log_file, quiet)
            download_singleton(url, output_path, log_file, quiet)

    return 0


def check_nr_md5s(data_dir, log_file, quiet):
    """Check integrity of all files in a dir with their paired .md5 files"""
    md5_files = list([p.resolve() for p in data_dir.glob("*.md5")])
    for md5_file in md5_files:
        data_file = md5_file.with_suffix("")
        check.check_md5_gz(data_file, md5_file, log_file, quiet)
    return 0


def process_nr(output_dir, date_prefix, log_file, quiet):

    nr_download(output_dir, date_prefix, log_file, quiet)
    check_nr_md5s(output_dir, log_file, quiet)

    tax_tar = "{}.{}".format(date_prefix, "taxdump.tar.gz")
    tax_tar_path = output_dir / pathlib.Path(tax_tar)

    with tarfile.open(tax_tar_path, "r:gz") as tar:
        for dmp in ["names.dmp", "nodes.dmp"]:
            message = "Extracting {} from taxdump.tar.gz".format(dmp)
            shared.give_user_feedback(message, log_file, quiet)
            tar.extract(dmp, path=output_dir)

            # Timestamp the extracted dmp file
            fout = output_dir / pathlib.Path(dmp)
            timestamped_fname = "{}.{}".format(date_prefix, dmp)
            timestamped_fout = output_dir / pathlib.Path(timestamped_fname)
            fout.rename(timestamped_fout)


def run():
    args = parse_arguments()
    if args.db == "nr":
        process_nr(args.output_dir, args.date, args.log_file, args.quiet)
    elif args.db == "gtdb":
        print("Not implemented yet")
        sys.exit(1)

