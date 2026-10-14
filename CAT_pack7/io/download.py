#!/usr/bin/env python3

import hashlib
import logging
import os
import pathlib
import shutil
import sys
import tarfile
import urllib.request
from collections import namedtuple
from typing import Callable

from CAT_pack7.config.settings import DownloadSettings
from CAT_pack7.config.validation import get_date
from CAT_pack7.io.parsers import optionally_compressed_handle
from CAT_pack7.utils import check
from CAT_pack7.utils.logging import Report, Status

log = logging.getLogger("CAT_pack")
DownloadReport = Callable[[pathlib.Path, int, int | None], None]
DownloadQueueReport = Callable[[list[pathlib.Path]], None]


def download_singleton(target_url, local_path: pathlib.Path, report: DownloadReport) -> None:
    """Download a single file to the specified location."""
    def reporthook(blocknum: int, blocksize: int, totalsize: int):
        total = totalsize if totalsize >= 0 else None
        downloaded = blocknum * blocksize
        if total is not None:
            downloaded = min(downloaded, total)
        report(local_path, downloaded, total)

    try:
        if report is not None:
            report(local_path, 0, None)
        urllib.request.urlretrieve(target_url, local_path, reporthook=reporthook)
    except Exception:
        log.info("Failed downloading file: {0}.".format(target_url))
        raise
        
    return





def multi_download(settings: DownloadSettings, report: Report,
    download_report: DownloadReport, queue_report: DownloadQueueReport
):
    """Download all required nr files in the specified output dir."""
    #existing_files = list([p.resolve() for p in settings.out_dir.iterdir()])

    urls = settings.urls.urls
    prefix = get_date()
    output_paths = {url.split('/')[-1]: settings.out_dir / f"{prefix}.{url.split('/')[-1]}" for url in urls}
    if queue_report is not None:
        queue_report(output_paths.())
    if report is not None:
        report("Downloading", Status.RUNNING, 0, len(urls))

    for index, (url, output_path) in enumerate(zip(urls, output_paths.values()), start=1):
        url_leaf = url.split("/")[-1]
        # @Bastiaan want to keep this functionality in?
        # if output_path in existing_files:
        #     log.info(
        #             "Skipping download of file {0}. It already exists.".format(
        #                 output_path.name
        #                 )
        #             )
        log.info("Downloading {0}.".format(url_leaf))
            
        download_singleton(url, output_path, download_report)
        if report is not None:
            report("Downloading", Status.RUNNING, index, len(urls))

    return output_paths


def check_nr_md5s(settings: DownloadSettings):
    """Check integrity of all files in a dir with their paired .md5 files.
    """
    md5_files = list([p.resolve() for p in settings.out_dir.glob("*.md5")])
    
    for md5_file in md5_files:
        data_file = md5_file.with_suffix("")
        check.check_md5_gz(data_file, md5_file)
                        
    return


def process_nr(output_dir, log_file, quiet, prefix, cleanup):
    nr_urls = [
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5",
        ("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
            "accession2taxid/prot.accession2taxid.FULL.gz"),
        ("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
            "accession2taxid/prot.accession2taxid.FULL.gz.md5"),
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5",
        ]
    
    # Fetch files.
    multi_download(nr_urls, output_dir, log_file, quiet, prefix)

    # Check files.
    check_nr_md5s(output_dir, log_file, quiet)

    # Process files.
    tax_tar = "{0}.{1}".format(prefix, "taxdump.tar.gz")
    tax_tar_path = output_dir / pathlib.Path(tax_tar)

    with tarfile.open(tax_tar_path, "r:gz") as tar:
        for dmp in ["names.dmp", "nodes.dmp"]:
            message = "Extracting {0} from taxdump.tar.gz".format(dmp)
            
            tar.extract(dmp, path=output_dir)

            # Timestamp the extracted dmp file
            outf1 = output_dir / pathlib.Path(dmp)
            timestamped_fname = "{0}.{1}".format(prefix, dmp)
            timestamped_outf1 = output_dir / pathlib.Path(timestamped_fname)
            outf1.rename(timestamped_outf1)

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
    
    log.info(
            "\n-----------------\n\n"
            "Done!\n\n"
            "A CAT_pack database can be build with:\n\n"
            "CAT_pack prepare \\\n"
            "--db_fasta {0} \\\n"
            "--names {1} \\\n"
            "--nodes {2} \\\n"
            "--acc2tax {3} \\\n"
            "--db_dir path/to/prepare_output\n".format(
                nr_gz.resolve(),
                names_dmp.resolve(),
                nodes_dmp.resolve(),
                acc2taxid_gz.resolve(),
                )
            )
    
    return

# GTDB.

## GENERAL.
prefixes_to_rank_names = {
        "d__": "superkingdom",  # Using superkingdom for compatibility with NCBI.
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

## FUNCTIONS.
def get_gtdb_latest_version(parent_url):
    """Read the version number from the VERSION file.
    """
    version_url = "{0}/VERSION.txt".format(parent_url)

    with urllib.request.urlopen(version_url) as f:
        version_data = f.read().decode()

    version = "{0} ({3})".format(*version_data.split("\n"))

    return version


def load_gtdb_md5sums(md5sums_file):
    """Create a dictionary from the MD5SUMS file."""
    md5_dict = {}
    
    with open(md5sums_file, "r") as f1:
        for line in f1:
            fields = [f.strip() for f in line.split()]
            fname = pathlib.Path(fields[1]).name
            md5_dict[fname] = fields[0]
            
    return md5_dict


def check_gtdb_md5s(data_dir, md5_dict, log_file, quiet):
    for f in data_dir.glob("*.gz"):
        if f.name not in md5_dict:
            continue
            
        md5 = check.gz_md5(f)
        if md5_dict[f.name] != md5:
            log.info("MD5 of {0} does not check out.".format(f.resolve()))

            sys.exit(1)
        else:
            log.info("MD5 of {0} checks out.".format(f.resolve()))
            
    return


def concatenate_taxonomies_gz(ar_tsv_gz, bac_tsv_gz, output_tsv):
    with open(output_tsv, "w") as outf1:
        for f in [ar_tsv_gz, bac_tsv_gz]:
            with optionally_compressed_handle(f, "r") as f1:
                for line in f1:
                    outf1.write(line)

    return


def parent_child_pairs(lineage_string):
    """Turn a ";" separated lineage to a list of tuples for child parent.

    Given a lineage string of the form

        "d__XXX;p__XXX;o__XXX;c__XXX;f__XXX;g_XXX;s__XXX"

    prepend "root" and split it into parent-child pairs

    [
        ("root", "d__XXX"),
        ("d__XXX", "p__XXX"),
        ("p__XXX", "o__XXX"),
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
    """Write the nodes.dmp from the taxonomy files."""
    seen_taxids = set()
    
    with open(taxonomies_tsv, "r") as f1, open(nodes_dmp, "w") as outf1:
        for line in f1:
            fields = [f.strip() for f in line.split("\t")]
            lineage_string = fields[1]
            pairs = parent_child_pairs(lineage_string)
            
            for pair in pairs:
                parent, child = pair[0], pair[1]
                
                if parent not in seen_taxids:
                    if parent == "root":
                        outf1.write(
                                "{0}{1}".format(
                                    "\t|\t".join(["root", "root", "no rank"]),
                                    "\t|\n"
                                    )
                                )
                        seen_taxids.add("root")
                    else:
                        seen_taxids.add(parent)
                        
                if child not in seen_taxids:
                    seen_taxids.add(child)
                    outf1.write(
                            "{0}{1}".format(
                                "\t|\t".join(
                                    [
                                        child,
                                        parent,
                                        prefixes_to_rank_names[child[:3]]
                                        ]
                                    ),
                                "\t|\n"
                                )
                            )
                    
    return


def write_names_dmp(taxonomies_tsv, names_dmp):
    """Write the names.dmp from the taxonomy files."""
    seen_taxids = set()
    
    with open(taxonomies_tsv, "r") as f1, open(names_dmp, "w") as outf1:
        outf1.write(
                "{0}{1}".format(
                    "\t|\t".join(["root", "root", "", "scientific name"]),
                    "\t|\n"
                    )
                )

        for line in f1:
            fields = [f.strip() for f in line.split("\t")]
            lineage_string = fields[1]

            for taxid in lineage_string.split(";"):
                if taxid not in seen_taxids:
                    outf1.write(
                            "{0}{1}".format(
                                "\t|\t".join(
                                    [
                                        taxid,
                                        taxid.split("__", 1)[-1],
                                        "",
                                        "scientific name"
                                        ]
                                    ),
                                "\t|\n"
                                )
                            )
                    seen_taxids.add(taxid)

    return


def genome_id_to_taxid(taxonomy_tsv):
    """Return a dictionary with GTDB taxid for each genome.
    """
    mapping = {}
    
    with open(taxonomy_tsv, "r") as f1:
        for line in f1:
            fields = [f.strip() for f in line.split("\t")]
            genome_id = fields[0]
            taxid = fields[1].split(";")[-1]
            mapping[genome_id] = taxid
            
    return mapping


def fastaIterator_gz(fasta_in_gz, gid2taxid):
    """Yield fastaRecord tuples with more information.

    This is adapted from biopython's SimpleFastaParser.
    https://github.com/biopython/biopython/blob/eec86d4bcb04bfcf86495f92f12faf3ff98a288d/Bio/SeqIO/FastaIO.py#L24

    The gid2taxid dictionary holds the mapping of a given genome accession
    (RS_GCF_XXXX or GB_GCA_XXXX). The taxid is propagated to all protein
    sequences of that genome.

    Positional argunents:
      fasta_in: pathlib.Path object: Path to the fasta.
      gid2taxid: dict: A dictionary with taxid for a genome accession of the
        form {"RS_CCF_XXXXX" : "s__Escherichia coli", ...}.

    Return:
      None, if the file is not a valid fasta. Breaks the iteration.

    Yields fastaRecord objects which are named tuples holding the following
    information:
      - id: str: Unique id of the fasta header, anything between the ">" and
        the first space.
      - seq: str: The sequence, capitalized and with trailing "*" stripped off
      - uid: str: Unique id, "_" joined md5sum and length.
      - taxid: str: GTDB taxonomy that was assigned to the genome.
    """
    origin = fasta_in_gz.name.replace("_protein.faa.gz", "")
    taxid = gid2taxid[origin]
    with optionally_compressed_handle(fasta_in_gz, "r") as f1:
        for line in f1:
            if line[0] == ">":
                title = line[1:].rstrip()
                break
            else:
                # No break encountered - probably an empty file.
                return
            
        lines = []
        for line in f1:
            if line[0] == ">":
                name = title.split(" ")[0]
                seq = "".join(
                        lines).replace(" ", "").replace("\r", "").rstrip("*")
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
    """Get a dictionary of duplicate uids (md5sum + length)."""
    seen_uids = {}
    multiplets = {}
    seq_counter, file_counter = 0, 0
    
    with optionally_compressed_handle(acc2taxid_fp, "w") as outf1:
        outf1.write("accession.version\ttaxid\n")
        
        for f in proteins_dir.rglob("*/*.faa.gz"):
            file_counter += 1
            
            for record in fastaIterator_gz(f, gid2taxid):
                # Write an entry.
                outf1.write("{0}\t{1}\n".format(record.id, record.taxid))

                # Check if duplicate.
                is_multiplet = False
                if record.uid not in seen_uids:
                    # Store the first occurrence id.
                    seen_uids[record.uid] = record.id
                else:
                    is_multiplet = True

                if is_multiplet is True:
                    if record.uid not in multiplets:
                        multiplets[record.uid] = fastaRecord(
                                "\x01".join(
                                    [seen_uids[record.uid], record.id]),
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
                log.info("Parsed {0:,d} sequences from {1:,d} files.".format(
                        seq_counter, file_counter))
        # This else is part of the outter for-loop.
        # It executes when the for loop finishes.
        else:
            # Create some whitespace for aligned printing.
            padding = len("[YYYY-MM-DD HH:MM:SS] ") * " "
            
            # Calculate the total number of identified multiplets.
            redundants = sum(map(len, [v for v in multiplets.values()]))
            
            log.info(
                    "    Total files: {0:>12,d}\n"
                    "{1}Total sequences: {2:>12,d}\n"
                    "{3}     Multiplets: {4:>12,d}\n"
                    "{5}of which unique: {6:>12,d}"
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
            
    return multiplets


def write_singletons(
    proteins_dir, duplicates, gid2taxid, singletons_fp, log_file, quiet):
    seq_counter, file_counter, skipped = 0, 0, 0
    
    with optionally_compressed_handle(singletons_fp, "w") as outf1:
        for f in proteins_dir.rglob("*/*.faa.gz"):
            file_counter += 1
            
            for record in fastaIterator_gz(f, gid2taxid):
                if record.uid not in duplicates:
                    outf1.write(">{0}\n{1}\n".format(record.id, record.seq))
                    seq_counter += 1
                else:
                    skipped += 1
                    
            if file_counter % 1000 == 0 and file_counter != 0:
                log.info("Written {0:,d} sequences from {1:,d} files "
                        "({2:,d} skipped).".format(
                            seq_counter, file_counter, skipped)
                        )
        else:
            log.info("Written {0:,d} sequences from {1:,d} files "
                    "({2:,d} skipped).".format(
                        seq_counter, file_counter, skipped)
                    )
            
    return


def concatenate_trees(bac_tree_fp, ar_tree_fp, all_tree_fp):
    """Concatenate the newick trees under a common root in a new file."""

    # Load the Bacteria tree as a string and make it a subtree.
    bac_tree = bac_tree_fp.read_text()
    bac_tree = bac_tree.rstrip().replace(
        "d__Bacteria;", "'100.0:d__Bacteria':1.0"
    )
    
    # Load the Archaea tree as a string and make it a subtree.
    ar_tree = ar_tree_fp.read_text()
    ar_tree = ar_tree.rstrip().replace("d__Archaea;", "'100.0:d__Archaea':1.0")
    
    # Concatenate the subtrees under a node named root.
    all_tree = "({0},{1})root;\n".format(ar_tree, bac_tree)
    
    # Write the file.
    all_tree_fp.write_text(all_tree)
    
    return


def process_gtdb(output_dir, log_file, quiet, cleanup=False):
    # Using "latest" as an entry point.
    # This needs to be checked for future versions.
    # Use the European mirror.
    parent_url = "https://data.gtdb.aau.ecogenomic.org/releases/latest/"

    version = get_gtdb_latest_version(parent_url)
    
    log.info("CAT_pack will download files from GTDB {0}.".format(version))
    

    file_list = [
            "VERSION.txt",
            "ar53_taxonomy.tsv.gz",
            "bac120_taxonomy.tsv.gz",
            "MD5SUM.txt",
            "bac120.tree",
            "ar53.tree",
            "genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"
            ]
    gtdb_urls = ["{0}/{1}".format(parent_url, file_) for file_ in file_list]

    # Fetch files.
    multi_download(gtdb_urls, output_dir, log_file, quiet, prefix=None)

    # Check files.
    md5sums_file = output_dir / pathlib.Path("MD5SUM.txt")
    md5sums_dict = load_gtdb_md5sums(md5sums_file)
    check_gtdb_md5s(output_dir, md5sums_dict, log_file, quiet)

    # Concatenate taxonomies.
    bacteria_tsv_gz = list(output_dir.glob("*bac*_taxonomy*"))[0]
    archaea_tsv_gz = list(output_dir.glob("*ar*_taxonomy*"))[0]
    all_taxa_tsv = output_dir / pathlib.Path("all_taxonomies.tsv")
    if not all_taxa_tsv.exists():
        concatenate_taxonomies_gz(
                archaea_tsv_gz, bacteria_tsv_gz, all_taxa_tsv)

    # Concatenate newick trees.
    bac_tree_fp = list(output_dir.glob("*bac*.tree"))[0]
    ar_tree_fp = list(output_dir.glob("*ar*.tree"))[0]
    concatenated_tree_fp = output_dir / pathlib.Path("gtdb.tree")
    if not concatenated_tree_fp.exists():
        log.info("Concatenating newick trees.")
        
        concatenate_trees(bac_tree_fp, ar_tree_fp, concatenated_tree_fp)

    # Extract protein files from archive.
    proteins_tar = list(output_dir.glob("*proteins_aa_reps*"))[0]
    proteins_dir = output_dir / pathlib.Path("protein_faa_reps")
    if not proteins_dir.is_dir():
        with tarfile.open(proteins_tar, "r:gz") as tar:
            # Fix for CVE-2007-4559.
            def is_within_directory(directory, target):
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(
                    tar, path=".", members=None, *, numeric_owner=False):
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner)
                
            safe_extract(tar, output_dir)
    else:
        log.info("Proteins directory {0} already exists.".format(
                proteins_dir.resolve()))
        

    # Process files.
    # NODES.
    nodes_dmp = output_dir / pathlib.Path("nodes.dmp")
    if not nodes_dmp.exists():
        log.info("Writing nodes information to {0}.".format(
                nodes_dmp.resolve()))
        
        write_nodes_dmp(all_taxa_tsv, nodes_dmp)
    else:
        log.info("Nodes file found : {0}.".format(nodes_dmp.resolve()))
        

    # NAMES.
    names_dmp = output_dir / pathlib.Path("names.dmp")
    if not names_dmp.exists():
        log.info("Writing names information to {0}.".format(
                names_dmp.resolve()))
        
        
        write_names_dmp(all_taxa_tsv, names_dmp)
    else:
        log.info("Names file found : {0}.".format(names_dmp.resolve()))
        

    gid2taxid = genome_id_to_taxid(all_taxa_tsv)

    # SEQUENCES.
    duplicates_fp = output_dir / pathlib.Path("dups.fa.gz")
    singletons_fp = output_dir / pathlib.Path("singletons.fa.gz")
    all_seqs_fp = output_dir / pathlib.Path("gtdb_seqs.fa.gz")
    acc2taxid_fp = output_dir / pathlib.Path("prot.accession2taxid.txt.gz")
    
    if not acc2taxid_fp.exists():
        log.info("1st pass: Extracting multiplets.")
        
        
        ## 1st pass.
        ## Write acc2taxid file, extract duplicates in their own fasta.
        duplicates = extract_duplicates(
            proteins_dir, gid2taxid, acc2taxid_fp, log_file, quiet)
        with optionally_compressed_handle(duplicates_fp, "w") as outf1:
            for rec in duplicates.values():
                outf1.write(">{0}\n{1}\n".format(rec.id, rec.seq))

        # 2nd pass.
        # Write the unique sequences to a separate file.
        log.info("2nd pass: Retrieving unique sequences.")
        
        
        write_singletons(
                proteins_dir,
                duplicates,
                gid2taxid,
                singletons_fp,
                log_file,
                quiet
                )

        log.info("Concatenating sequence files.")
        
        
        # Concatenate the two files into one.
        with optionally_compressed_handle(all_seqs_fp, "w") as outf1:
            for f in [duplicates_fp, singletons_fp]:
                with optionally_compressed_handle(f, "r") as f1:
                    for line in f1:
                        outf1.write(line)

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
        log.info("Cleanup specified. Removing unnecessary files and folders.")
        

        for target in remove_targets:
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
                
    log.info(
            "\n-----------------\n\n"
            "Done!\n\n"
            "A CAT_pack database can be build with:\n\n"
            "CAT_pack prepare \\\n"
            "--db_fasta {0} \\\n"
            "--names {1} \\\n"
            "--nodes {2} \\\n"
            "--acc2tax {3} \\\n"
            "--db_dir path/to/prepare_output\n".format(
                all_seqs_fp.resolve(),
                names_dmp.resolve(),
                nodes_dmp.resolve(),
                acc2taxid_fp.resolve(),
                )
            )


    
    return

if __name__ == "__main__":
    sys.exit("Run \'CAT_pack download\' to download and preprocess data from "
            "NCBI nr or GTDB.")
