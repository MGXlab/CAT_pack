import argparse
from collections import namedtuple
import datetime
import hashlib
import pathlib
import shutil
import sys
import tarfile
import urllib.request
import urllib.parse

import shared
import check


def parse_arguments():
    date = datetime.datetime.now().strftime("%Y-%m-%d")

    parser = argparse.ArgumentParser(
        prog="CAT download",
        description=(
            "Download and preprocess sequence and taxonomy information. "
            "Currently supports the NCBI nr and GTDB"
        ),
        usage="CAT download -db [nr|gtdb] -o output_dir",
        add_help=False,
    )

    required = parser.add_argument_group("Required Arguments")
    shared.add_argument(required, "db", True)
    shared.add_argument(required, "output_dir", True)

    optional = parser.add_argument_group("Optional Arguments")
    shared.add_argument(optional, "cleanup", False)
    shared.add_argument(optional, "quiet", False)
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


def multi_download(url_list, output_dir, log_file, quiet, prefix=None):
    """Download all required nr files in the specified output dir"""

    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        existing_files = []
    else:
        existing_files = list([p.resolve() for p in output_dir.iterdir()])

    for url in url_list:
        url_leaf = url.split("/")[-1]
        if prefix:
            output_basename = "{}.{}".format(prefix, url_leaf)
        else:
            output_basename = url_leaf
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


def process_nr(output_dir, log_file, quiet, prefix, cleanup):
    nr_urls = [
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz.md5",
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5",
    ]
    # Fetch files
    multi_download(nr_urls, output_dir, log_file, quiet, prefix)

    # Check files
    check_nr_md5s(output_dir, log_file, quiet)

    # Process files
    tax_tar = "{}.{}".format(prefix, "taxdump.tar.gz")
    tax_tar_path = output_dir / pathlib.Path(tax_tar)

    with tarfile.open(tax_tar_path, "r:gz") as tar:
        for dmp in ["names.dmp", "nodes.dmp"]:
            message = "Extracting {} from taxdump.tar.gz".format(dmp)
            shared.give_user_feedback(message, log_file, quiet)
            tar.extract(dmp, path=output_dir)

            # Timestamp the extracted dmp file
            fout = output_dir / pathlib.Path(dmp)
            timestamped_fname = "{}.{}".format(prefix, dmp)
            timestamped_fout = output_dir / pathlib.Path(timestamped_fname)
            fout.rename(timestamped_fout)

    if cleanup is True:
        targets = [tax_tar_path.resolve()]
        for i in output_dir.glob("*.md5"):
            targets.append(i.resolve())
        for t in targets:
            t.unlink()

    # CAT prepare
    nr_gz = list(output_dir.glob("*nr.gz"))[0]
    names_dmp = list(output_dir.glob("*names.dmp"))[0]
    nodes_dmp = list(output_dir.glob("*nodes.dmp"))[0]
    acc2taxid_gz = list(output_dir.glob("*accession2taxid.FULL.gz"))[0]
    message = (
        "\n-----------------\n\n"
        "Done!\n"
        "A CAT database can be build with\n\n"
        "CAT prepare \\\n"
        "--db_fasta {} \\\n"
        "--names {} \\\n"
        "--nodes {} \\\n"
        "--acc2tax {} \\\n"
        "--db_dir path/to/prepare_output\n".format(
            nr_gz.resolve(),
            names_dmp.resolve(),
            nodes_dmp.resolve(),
            acc2taxid_gz.resolve(),
        )
    )
    shared.give_user_feedback(message, log_file, quiet, show_time=False)


# GTDB

## GENERAL
prefixes_to_rank_names = {
    "d__": "superkingdom",  # Using superkingdom for compatibility with ncbi
    "p__": "phylum",
    "o__": "order",
    "c__": "class",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}

fastaRecord = namedtuple(
    "fastaRecord",
    ["id", "seq", "uid", "taxid"],
)

## FUNCTIONS
def get_gtdb_latest_version():
    """Read the version number from the VERSION file"""
    version_url = "https://data.gtdb.ecogenomic.org/releases/latest/VERSION"

    with urllib.request.urlopen(version_url) as f:
        version_data = f.read().decode()

    version = version_data.split("\n")[0].lstrip("v")

    return version


def load_gtdb_md5sums(md5sums_file):
    """Create a dictionary from the MD5SUMS file"""
    md5_dict = {}
    with open(md5sums_file, "r") as fin:
        for line in fin:
            fields = [f.strip() for f in line.split()]
            fname = pathlib.Path(fields[1]).name
            md5_dict[fname] = fields[0]

        # at the time of GTDB release v207, and looking at earlier versions,
        # the filenames in the "latest" section (https://data.gtdb.ecogenomic.org/releases/latest/)
        # look like this: "bac120_taxonomy.tsv.gz", but the names
        # in the MD5SUM file from the same location look like 
        # this: "bac120_taxonomy_r207.tsv.gz"
        # So, here doing an ad hoc check and removing
            # getting listed version
            # seeing if that's in the file name as _r<version>
            # and removing that string from our dict keys if so

        # getting version
        gtdb_version = get_gtdb_latest_version()
        version_string = f"_r{gtdb_version}"

        for key in md5_dict:

            if version_string in key:

                # removing version string
                new_name = key.replace(version_string, "")

                md5_dict[new_name] = md5_dict.pop(key)

    return md5_dict


def check_gtdb_md5s(data_dir, md5_dict, log_file, quiet):
    for f in data_dir.glob("*.gz"):
        md5 = check.gz_md5(f)
        if md5_dict[f.name] != md5:
            message = "MD5 of {0} does not check out.".format(f.resolve())
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
        else:
            message = "MD5 of {0} checks out.".format(f.resolve())
            shared.give_user_feedback(message, log_file, quiet)
    return


def concatenate_taxonomies_gz(ar_tsv_gz, bac_tsv_gz, output_tsv):
    with open(output_tsv, "w") as fout:
        for f in [ar_tsv_gz, bac_tsv_gz]:
            with shared.optionally_compressed_handle(f, "r") as fin:
                for line in fin:
                    fout.write(line)

    return


def parent_child_pairs(lineage_string):
    """Turn a ";" separated lineage to a list of tuples for child parent

    Given a lineage string of the form

        'd__XXX;p__XXX;o__XXX;c__XXX;f__XXX;g_XXX;s__XXX'

    prepend 'root' and split it into parent-child pairs

    [
        ('root', 'd__XXX'),
        ('d__XXX', 'p__XXX'),
        ('p__XXX', 'o__XXX'),
        ...
    ]
    """
    with_root = ";".join(["root", lineage_string])
    lineage_list = with_root.split(";")
    pairs = []
    for i in range(0, len(lineage_list) - 1):
        parent, child = lineage_list[i], lineage_list[i + 1]
        pairs.append((parent, child))
    return pairs


def write_nodes_dmp(taxonomies_tsv, nodes_dmp):
    """Write the nodes.dmp from the taxonomies files"""
    seen_taxids = []
    with open(taxonomies_tsv, "r") as fin, open(nodes_dmp, "w") as fout:
        for line in fin:
            fields = [f.strip() for f in line.split("\t")]
            lineage_string = fields[1]
            pairs = parent_child_pairs(lineage_string)
            for pair in pairs:
                parent, child = pair[0], pair[1]
                if parent not in seen_taxids:
                    if parent == "root":
                        fout.write(
                            "{}{}".format(
                                "\t|\t".join(["root", "root", "no rank"]),
                                "\t|\n",
                            )
                        )
                        seen_taxids.append("root")
                    else:
                        seen_taxids.append(parent)
                if child not in seen_taxids:
                    seen_taxids.append(child)
                    fout.write(
                        "{}{}".format(
                            "\t|\t".join(
                                [
                                    child,
                                    parent,
                                    prefixes_to_rank_names[child[:3]],
                                ]
                            ),
                            "\t|\n",
                        )
                    )
    return


def write_names_dmp(taxonomies_tsv, names_dmp):
    seen_taxids = []
    with open(taxonomies_tsv, "r") as fin, open(names_dmp, "w") as fout:
        fout.write(
            "{}{}".format(
                "\t|\t".join(["root", "root", "scientific name"]), "\t|\n"
            )
        )
        for line in fin:
            taxid = line.split(";")[-1].strip()
            if taxid not in seen_taxids:
                fout.write(
                    "{}{}".format(
                        "\t|\t".join([taxid, taxid, "scientific name"]),
                        "\t|\n",
                    )
                )
                seen_taxids.append(taxid)
    return


def genome_id_to_taxid(taxonomy_tsv):
    """Return a dictionary with GTDB taxid for each genome"""
    mapping = {}
    with open(taxonomy_tsv, "r") as fin:
        for line in fin:
            fields = [f.strip() for f in line.split("\t")]
            genome_id = fields[0]
            taxid = fields[1].split(";")[-1]
            mapping[genome_id] = taxid
    return mapping


def fastaIterator(fasta_in, gid2taxid):
    """Yield fastaRecord tuples with more information

    This is adapted from biopython's SimpleFastaParser.
    https://github.com/biopython/biopython/blob/eec86d4bcb04bfcf86495f92f12faf3ff98a288d/Bio/SeqIO/FastaIO.py#L24

    The gid2taxid dictionary holds the mapping of a given genome accession
    (RS_GCF_XXXX or GB_GCA_XXXX). The taxid is propagated to all protein
    sequences of that genome.


    Positional argunents:
      fasta_in: pathlib.Path object: Path to the fasta
      gid2taxid: dict: A dictionart with taxid for a genome accession of the
        form {"RS_CCF_XXXXX" : "s__Escherichia coli", ...}

    Return:
      None, if the file is not a valid fasta. Breaks the iteration

    Yields fastaRecord objects which are named tuples holding the following
    information:
      - id: str: Unique id of the fasta header, anything between the '>' and the
        first space
      - seq: str: The sequence, capitalized and with trailing '*' stripped off
      - uid: str: Unique id, '_' joined md5sum and length
      - taxid: str: GTDB taxonomy that was assigned to the genome
    """
    origin = fasta_in.name.replace("_protein.faa", "")
    taxid = gid2taxid[origin]
    with open(fasta_in, "r") as fin:
        for line in fin:
            if line[0] == ">":
                title = line[1:].rstrip()
                break
            else:
                # no break encountered - probably an empty file
                return
        lines = []
        for line in fin:
            if line[0] == ">":
                name = title.split(" ")[0]
                seq = (
                    "".join(lines)
                    .replace(" ", "")
                    .replace("\r", "")
                    .rstrip("*")
                )
                length = len(seq)
                md5sum = hashlib.md5(seq.encode()).hexdigest()
                uid = "_".join([md5sum, str(length)])
                yield fastaRecord(name, seq, uid, taxid)
                lines = []
                title = line[1:].rstrip()
                continue
            lines.append(line.rstrip())

        name = title.split(" ")[0]
        seq = "".join(lines).replace(" ", "").replace("\r", "").rstrip("*")
        length = len(seq)
        md5sum = hashlib.md5(seq.encode()).hexdigest()
        uid = "_".join([md5sum, str(length)])
        yield fastaRecord(name, seq, uid, taxid)


def extract_duplicates(proteins_dir, gid2taxid, acc2taxid_fp, log_file, quiet):
    """Get a dictionary of duplicate uids (md5sum + length)"""
    seen_uids = {}
    multiplets = {}
    seq_counter, file_counter = 0, 0
    with shared.optionally_compressed_handle(acc2taxid_fp, "w") as fout:
        fout.write("accession.version\ttaxid\n")
        for f in proteins_dir.rglob("*/*.faa"):
            file_counter += 1
            for record in fastaIterator(f, gid2taxid):
                # Write an entry
                fout.write("{}\t{}\n".format(record.id, record.taxid))

                # Check if duplicate
                is_multiplet = False
                if record.uid not in seen_uids:
                    # Store the first occurrence id
                    seen_uids[record.uid] = record.id
                else:
                    is_multiplet = True

                if is_multiplet is True:
                    if record.uid not in multiplets:
                        multiplets[record.uid] = fastaRecord(
                            "\x01".join([seen_uids[record.uid], record.id]),
                            record.seq,
                            None,
                            None,
                        )
                    else:
                        old_rec = multiplets[record.uid]
                        new_rec = fastaRecord(
                            "\x01".join([old_rec.id, record.id]),
                            old_rec.seq,
                            None,
                            None,
                        )

                        multiplets[record.uid] = new_rec

                seq_counter += 1

            if file_counter % 1000 == 0 and file_counter != 0:
                message = "Parsed {} sequences from {} files".format(
                    seq_counter, file_counter
                )
                shared.give_user_feedback(message, log_file, quiet)
        # This else is part of the outter for-loop
        # It executes when the for loop finishes
        else:
            # Create some whitespace for aligned printing
            padding = len("[YYYY-MM-DD HH:MM:SS] ") * " "
            # Calculate the total number of identified multiplets
            redundants = sum(map(len, [v for v in multiplets.values()]))
            message = (
                "    Total files: {:>12}\n"
                "{}Total sequences: {:>12}\n"
                "{}     Multiplets: {:>12}\n"
                "{}of which unique: {:>12}"
                "".format(
                    file_counter,
                    padding,
                    seq_counter,
                    padding,
                    redundants,
                    padding,
                    len(multiplets),
                )
            )
            shared.give_user_feedback(message, log_file, quiet)
    return multiplets


def write_singletons(
    proteins_dir, duplicates, gid2taxid, singletons_fp, log_file, quiet
):
    seq_counter, file_counter, skipped = 0, 0, 0
    with shared.optionally_compressed_handle(singletons_fp, "w") as fout:
        for f in proteins_dir.rglob("*/*.faa"):
            file_counter += 1
            for record in fastaIterator(f, gid2taxid):
                if record.uid not in duplicates:
                    fout.write(">{}\n{}\n".format(record.id, record.seq))
                    seq_counter += 1
                else:
                    skipped += 1
            if file_counter % 1000 == 0 and file_counter != 0:
                message = (
                    "Written {} sequences from {} files ({} skipped)".format(
                        seq_counter, file_counter, skipped
                    )
                )
                shared.give_user_feedback(message, log_file, quiet)
        else:
            message = "Written {} sequences from {} files ({} skipped)".format(
                seq_counter, file_counter, skipped
            )
            shared.give_user_feedback(message, log_file, quiet)


def concatenate_trees(bac_tree_fp, ar_tree_fp, all_tree_fp):
    """Concatenate the newick trees under a common root in a new file"""

    # Load the bacteria tree as a string and make it a subtree
    bac_tree = bac_tree_fp.read_text()
    bac_tree = bac_tree.rstrip().replace(
        "d__Bacteria;", "'100.0:d__Bacteria':1.0"
    )
    # Load the Archaea tree as a string and make it a subtree
    ar_tree = ar_tree_fp.read_text()
    ar_tree = ar_tree.rstrip().replace("d__Archaea;", "'100.0:d__Archaea':1.0")
    # Concatenate the subtrees under a node named root
    all_tree = "({},{})root;\n".format(ar_tree, bac_tree)
    # Write the file
    all_tree_fp.write_text(all_tree)
    return


def process_gtdb(output_dir, log_file, quiet, cleanup=False):

    # Using `latest` as an entry point
    # This needs to be checked for future versions
    gtdb_urls = [
        "https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv.gz",
        "https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz",
        "https://data.gtdb.ecogenomic.org/releases/latest/MD5SUM",
        "https://data.gtdb.ecogenomic.org/releases/latest/bac120.tree",
        "https://data.gtdb.ecogenomic.org/releases/latest/ar53.tree",
        "https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz",
    ]

    # Fetch files
    multi_download(gtdb_urls, output_dir, log_file, quiet, prefix=None)

    # Check files
    md5sums_file = output_dir / pathlib.Path("MD5SUM")
    md5sums_dict = load_gtdb_md5sums(md5sums_file)

    check_gtdb_md5s(output_dir, md5sums_dict, log_file, quiet)

    # Concatenate taxonomies
    bacteria_tsv_gz = list(output_dir.glob("*bac*_taxonomy*"))[0]
    archaea_tsv_gz = list(output_dir.glob("*ar*_taxonomy*"))[0]
    all_taxa_tsv = output_dir / pathlib.Path("all_taxonomies.tsv")
    if not all_taxa_tsv.exists():
        concatenate_taxonomies_gz(
            archaea_tsv_gz, bacteria_tsv_gz, all_taxa_tsv
        )

    # Concatenate newick trees
    bac_tree_fp = list(output_dir.glob("*bac*.tree"))[0]
    ar_tree_fp = list(output_dir.glob("*ar*.tree"))[0]
    concatenated_tree_fp = output_dir / pathlib.Path("gtdb.tree")
    if not concatenated_tree_fp.exists():
        message = "Concatenating newick trees"
        shared.give_user_feedback(message, log_file, quiet)
        concatenate_trees(bac_tree_fp, ar_tree_fp, concatenated_tree_fp)

    # Extract protein files from archive
    proteins_tar = list(output_dir.glob("*proteins_aa_reps*"))[0]
    proteins_dir = output_dir / pathlib.Path("protein_faa_reps")
    if not proteins_dir.is_dir():
        with tarfile.open(proteins_tar, "r:gz") as tar:
            tar.extractall(output_dir)
    else:
        message = "Proteins directory {} already exists".format(
            proteins_dir.resolve()
        )
        shared.give_user_feedback(message, log_file, quiet)

    # Process files
    # NODES
    nodes_dmp = output_dir / pathlib.Path("nodes.dmp")
    if not nodes_dmp.exists():
        message = "Writing nodes information to {}".format(nodes_dmp.resolve())
        shared.give_user_feedback(message, log_file, quiet)
        write_nodes_dmp(all_taxa_tsv, nodes_dmp)
    else:
        message = "Nodes file found : {}".format(nodes_dmp.resolve())
        shared.give_user_feedback(message, log_file, quiet)

    # NAMES
    names_dmp = output_dir / pathlib.Path("names.dmp")
    if not names_dmp.exists():
        message = "Writing names information to {}".format(names_dmp.resolve())
        shared.give_user_feedback(message, log_file, quiet)
        write_names_dmp(all_taxa_tsv, names_dmp)
    else:
        message = "Names file found : {}".format(names_dmp.resolve())
        shared.give_user_feedback(message, log_file, quiet)

    gid2taxid = genome_id_to_taxid(all_taxa_tsv)

    # SEQUENCES
    duplicates_fp = output_dir / pathlib.Path("dups.fa.gz")
    singletons_fp = output_dir / pathlib.Path("singletons.fa.gz")
    all_seqs_fp = output_dir / pathlib.Path("gtdb_seqs.fa.gz")
    acc2taxid_fp = output_dir / pathlib.Path("prot.accession2taxid.txt.gz")
    if not acc2taxid_fp.exists():
        message = "1st pass: Extracting multiplets"
        shared.give_user_feedback(message, log_file, quiet)
        ## 1st pass
        ## Write acc2taxid file, extract duplicates in their own fasta
        duplicates = extract_duplicates(
            proteins_dir, gid2taxid, acc2taxid_fp, log_file, quiet
        )
        with shared.optionally_compressed_handle(duplicates_fp, "w") as fa:
            for rec in duplicates.values():
                fa.write(">{}\n{}\n".format(rec.id, rec.seq))

        # 2nd pass
        # Write the unique sequences to a separate file
        message = "2nd pass: Retrieving unique sequences"
        shared.give_user_feedback(message, log_file, quiet)
        write_singletons(
            proteins_dir, duplicates, gid2taxid, singletons_fp, log_file, quiet
        )

        message = "Concatenating sequence files"
        shared.give_user_feedback(message, log_file, quiet)
        # Concatenate the two files into one
        with shared.optionally_compressed_handle(all_seqs_fp, "w") as fout:
            for f in [duplicates_fp, singletons_fp]:
                with shared.optionally_compressed_handle(f, "r") as fin:
                    for line in fin:
                        fout.write(line)

    if cleanup is True:
        remove_targets = [
            proteins_dir,
            bac_tree_fp,
            ar_tree_fp,
            duplicates_fp,
            singletons_fp,
            bacteria_tsv_gz,
            archaea_tsv_gz,
            all_taxa_tsv,
            proteins_tar,
        ]
        message = "Cleanup specified. Removing unnecessary files and dirs"
        shared.give_user_feedback(message, log_file, quiet)

        for target in remove_targets:
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
    message = (
        "\n-----------------\n\n"
        "Done!\n"
        "A CAT database can be build with\n\n"
        "CAT prepare \\\n"
        "--db_fasta {} \\\n"
        "--names {} \\\n"
        "--nodes {} \\\n"
        "--acc2tax {} \\\n"
        "--db_dir path/to/prepare_output\n".format(
            all_seqs_fp.resolve(),
            names_dmp.resolve(),
            nodes_dmp.resolve(),
            acc2taxid_fp.resolve(),
        )
    )

    shared.give_user_feedback(message, log_file, quiet, show_time=False)


def run():
    args = parse_arguments()

    if args.no_log:
        log_file = None
    else:
        log_fname = "{}.CAT_download.log".format(args.date)
        log_file = args.output_dir / pathlib.Path(log_fname)

    setattr(args, "log_file", log_file)

    if args.db == "nr":
        process_nr(
            args.output_dir,
            args.log_file,
            args.quiet,
            prefix=args.date,
            cleanup=args.cleanup,
        )
    elif args.db == "gtdb":
        process_gtdb(args.output_dir, args.log_file, args.quiet, args.cleanup)
