#!/usr/bin/env python3

import gzip
import hashlib
import logging
import pathlib
import re
import shutil
import tarfile
import urllib.request
from collections import namedtuple
from contextlib import ExitStack, contextmanager, nullcontext
from functools import singledispatch
from itertools import chain
from typing import Callable

import io
from CAT_pack.config.settings import DownloadOutputs, DownloadSettings, GTDBDownloadOutputs
from CAT_pack.io.parsers import optionally_compressed_handle
from CAT_pack.utils import check
from CAT_pack.utils.errors import CheckSumError
from CAT_pack.utils.logging import Report, Status

log = logging.getLogger("CAT_pack")
DownloadReport = Callable[[pathlib.Path, int, int | None], None]
DownloadQueueReport = Callable[[list[pathlib.Path]], None]


def download_singleton(target_url: str, local_path: pathlib.Path, report: DownloadReport = lambda *_: None) -> None:
    """Download a single file to the specified location."""
    partial = local_path.with_name(f"{local_path.name}.part")

    def reporthook(blocknum: int, blocksize: int, totalsize: int):
        total = totalsize if totalsize >= 0 else None
        downloaded = blocknum * blocksize
        if total is not None:
            downloaded = min(downloaded, total)
        report(local_path, downloaded, total)

    try:
        partial.unlink(missing_ok=True)
        report(local_path, 0, None)
        urllib.request.urlretrieve(target_url, partial, reporthook=reporthook)
        size = partial.stat().st_size
        partial.replace(local_path)
        report(local_path, size, size)
    finally:
        partial.unlink(missing_ok=True)



def multi_download(settings: DownloadSettings, report: Report,
    download_report: DownloadReport, queue_report: DownloadQueueReport
):
    """Download all required nr or GTDB files in the specified output dir."""
    queue_report([item.path for item in settings.files.downloads])
    report("Downloading", Status.RUNNING, 0, len(settings.files.downloads))
    for index, item in enumerate(settings.files.downloads, 1):
        log.info("Downloading %s.", item.path.name)
        download_singleton(item.url, item.path, download_report)
        report("Downloading", Status.RUNNING, index, len(settings.files.downloads))

# FIXME: move to check.py
def check_nr_md5s(files: DownloadOutputs, report: Report) -> None:
    pairs = [item for item in files.downloads if item.checksum is not None]
    for index, item in enumerate(pairs, 1):
        check.check_md5_gz(item.path, item.checksum)
        report("Check files", Status.RUNNING, index, len(pairs))

# FIXME: move to parsers.py
def extract_nr_taxonomy(files: DownloadOutputs, report: Report) -> None:
    targets = (("names.dmp", files.names), ("nodes.dmp", files.nodes))
    with tarfile.open(files.taxonomy_archive, "r:gz") as archive:
        members = []
        for name, target in targets:

            member = archive.getmember(name)
            members.append((member, target))
        for index, (member, target) in enumerate(members, 1):
            log.info("Extracting %s to %s.", member.name, target)
            partial = target.with_name(f"{target.name}.part")
            try:
                partial.unlink(missing_ok=True)
                with archive.extractfile(member) as source, partial.open("wb") as output:
                    shutil.copyfileobj(source, output)
                partial.replace(target)
            finally:
                partial.unlink(missing_ok=True)
            report("Process files", Status.RUNNING, index, len(targets))


def cleanup_nr_downloads(files: DownloadOutputs) -> None:
    targets = [files.taxonomy_archive]
    targets.extend(item.checksum for item in files.downloads if item.checksum is not None)
    for target in targets:
        log.info(f"Removing {target}")
        target.unlink()

# FIXME: move to check.py
@singledispatch
def check_files(files: DownloadOutputs, report: Report) -> None:
    check_nr_md5s(files, report)

# FIXME: move to parsers.py?
@singledispatch
def process_files(files: DownloadOutputs, report: Report) -> None:
    extract_nr_taxonomy(files, report)


@singledispatch
def cleanup_downloads(files: DownloadOutputs) -> None:
    cleanup_nr_downloads(files)


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


@contextmanager
def staged_outputs(paths: list[pathlib.Path]):
    """Publish completed outputs only after every writer has closed successfully."""
    partials = [path.with_name(f"{path.name}.part") for path in paths]
    try:
        yield partials
        for partial, target in zip(partials, paths):
            partial.replace(target)
    finally:
        for partial in partials:
            partial.unlink(missing_ok=True)

# FIXME: move to parsers.py
def load_gtdb_md5sums(md5sums_file: pathlib.Path) -> dict[str, str]:
    """Create a dictionary from the MD5SUMS file."""
    md5_dict = {}
    for line in md5sums_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split(maxsplit=1)
        fname = pathlib.PurePosixPath(fields[1].lstrip("*")).name
        fname = re.sub(r"_r[0-9]+(?=\.)", "", fname)
        md5_dict[fname] = fields[0].lower()
    return md5_dict

# FIXME: move to check.py
@check_files.register
def check_gtdb_md5s(files: GTDBDownloadOutputs, report: Report) -> None:
    md5_dict = load_gtdb_md5sums(files.checksums)
    inputs = [item for item in files.downloads
              if item.path.suffix == ".gz" and item.url.rsplit("/", 1)[-1] in md5_dict]
    for index, item in enumerate(inputs, 1):
        fname = item.url.rsplit("/", 1)[-1]
        md5_exp = md5_dict[fname]
        md5 = check.gz_md5(item.path, block_size=1024 * 1024)
        if md5 != md5_exp:
            raise CheckSumError(f"MD5 of {fname} does not check out.", path=item.path)
        log.info("MD5 of %s checks out.", item.path.resolve())
        report("Check files", Status.RUNNING, index, len(inputs))
    version = " ".join(files.version.read_text(encoding="utf-8").split())
    log.info("Downloaded GTDB %s.", version)


@contextmanager
def _taxonomy_handle(taxonomies_tsv):
    if isinstance(taxonomies_tsv, (str, pathlib.Path)):
        taxonomies_tsv = (taxonomies_tsv,)
    with ExitStack() as stack:
        handles = [stack.enter_context(optionally_compressed_handle(path, "r"))
                   for path in taxonomies_tsv]
        yield chain.from_iterable(handles)


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

# FIXME: move to writers.py
def write_nodes_dmp(taxonomies_tsv, nodes_dmp):
    """Write the nodes.dmp from the taxonomy files."""
    seen_taxids = set()

    with _taxonomy_handle(taxonomies_tsv) as f1, open(nodes_dmp, "w") as outf1:
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

# FIXME: move to writers.py
def write_names_dmp(taxonomies_tsv, names_dmp):
    """Write the names.dmp from the taxonomy files."""
    seen_taxids = set()

    with _taxonomy_handle(taxonomies_tsv) as f1, open(names_dmp, "w") as outf1:
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

# FIXME: move to parsers.py
def genome_id_to_taxid(taxonomy_tsv):
    """Return a dictionary with GTDB taxid for each genome.
    """
    mapping = {}

    with _taxonomy_handle(taxonomy_tsv) as f1:
        for line in f1:
            fields = [f.strip() for f in line.split("\t")]
            genome_id = fields[0]
            taxid = fields[1].split(";")[-1]
            mapping[genome_id] = taxid

    return mapping

# FIXME: move to writers.py
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
    with staged_outputs([all_tree_fp]) as (partial,):
        partial.write_text(all_tree, encoding="utf-8")

    return


fastaRecord = namedtuple("fastaRecord", ["id", "seq", "uid", "taxid"])

# FIXME: move to parsers.py
# FIXME: remove duplicate code...
def fastaIterator_gz(fasta_in_gz, gid2taxid, *, handle=None):
    origin = pathlib.Path(fasta_in_gz).name.removesuffix(".gz").removesuffix("_protein.faa")
    taxid = gid2taxid[origin]
    source = nullcontext(handle) if handle is not None else optionally_compressed_handle(fasta_in_gz, "r")
    with source as f1:
        title = None
        lines = []
        for line in f1:
            if not line.strip():
                continue
            if line.startswith(">"):
                if title is not None:
                    name = title.split(" ")[0]
                    seq = "".join(lines).upper().rstrip("*")
                    uid = (hashlib.md5(seq.encode()).digest(), len(seq))
                    yield fastaRecord(f"{origin}~{name}", seq, uid, taxid)
                title = line[1:].rstrip()
                lines = []
            elif title is None:
                return
            else:
                lines.append("".join(line.split()))
        if title is not None:
            name = title.split(" ")[0]
            seq = "".join(lines).upper().rstrip("*")
            uid = (hashlib.md5(seq.encode()).digest(), len(seq))
            yield fastaRecord(f"{origin}~{name}", seq, uid, taxid)

# FIXME: move to parsers.py
def iter_gtdb_proteins(proteins_tar: pathlib.Path, gid2taxid: dict[str, str]):
    with tarfile.open(proteins_tar, "r|gz") as archive:
        for member in archive:
            name = pathlib.PurePosixPath(member.name).name
            if not name.endswith(("_protein.faa.gz", "_protein.faa")):
                continue
            with ExitStack() as stack:
                raw = stack.enter_context(archive.extractfile(member))
                source = stack.enter_context(gzip.GzipFile(fileobj=raw)) if name.endswith(".gz") else raw
                handle = stack.enter_context(io.TextIOWrapper(source, encoding="utf-8"))
                yield from fastaIterator_gz(pathlib.Path(name), gid2taxid, handle=handle)

# FIXME: move to writers.py and or parsers.py?
def extract_duplicates(proteins_tar, gid2taxid, acc2taxid_fp):
    """Write the accession table and return duplicate sequence records."""
    seen_uids = {}
    multiplets = {}
    seq_counter = 0
    with gzip.open(acc2taxid_fp, "wt", encoding="utf-8") as outf1:
        outf1.write("accession.version\ttaxid\n")
        for record in iter_gtdb_proteins(proteins_tar, gid2taxid):
            outf1.write(f"{record.id}\t{record.taxid}\n")
            if record.uid not in seen_uids:
                seen_uids[record.uid] = record.id
            elif record.uid not in multiplets:
                multiplets[record.uid] = fastaRecord(
                    "\x01".join([seen_uids[record.uid], record.id]), record.seq, None, None,
                )
            else:
                old_rec = multiplets[record.uid]
                multiplets[record.uid] = fastaRecord(
                    "\x01".join([old_rec.id, record.id]), old_rec.seq, None, None,
                )
            seq_counter += 1
            if seq_counter % 100000 == 0:
                log.info(f"Indexed {seq_counter:,} GTDB proteins.")
    log.info(f"Indexed {seq_counter:,} GTDB proteins in total.")
    return multiplets

# FIXME: move to writers.py
def write_singletons(proteins_tar, duplicates, gid2taxid, singletons_fp, *, append=False):
    with gzip.open(singletons_fp, "at" if append else "wt", encoding="utf-8") as outf1:
        for record in iter_gtdb_proteins(proteins_tar, gid2taxid):
            if record.uid not in duplicates:
                outf1.write(f">{record.id}\n{record.seq}\n")


@process_files.register
def process_gtdb(files: GTDBDownloadOutputs, report: Report) -> None:
    taxonomies_tsv = files.taxonomy_tables
    with staged_outputs([files.nodes, files.names]) as (nodes_dmp, names_dmp):
        write_nodes_dmp(taxonomies_tsv, nodes_dmp)
        write_names_dmp(taxonomies_tsv, names_dmp)
    gid2taxid = genome_id_to_taxid(taxonomies_tsv)
    report("Process files", Status.RUNNING, 1, 4)
    ar_tree_fp, bac_tree_fp = files.trees
    concatenate_trees(bac_tree_fp, ar_tree_fp, files.tree)
    report("Process files", Status.RUNNING, 2, 4)
    proteins_tar = files.protein_archive
    with staged_outputs([files.db_fasta, files.acc2tax]) as (all_seqs_fp, acc2taxid_fp):
        duplicates = extract_duplicates(proteins_tar, gid2taxid, acc2taxid_fp)
        report("Process files", Status.RUNNING, 3, 4)
        with gzip.open(all_seqs_fp, "wt", encoding="utf-8") as outf1:
            for rec in duplicates.values():
                outf1.write(f">{rec.id}\n{rec.seq}\n")
        write_singletons(proteins_tar, duplicates, gid2taxid, all_seqs_fp, append=True)
    report("Process files", Status.RUNNING, 4, 4)


@cleanup_downloads.register
def cleanup_gtdb_downloads(files: GTDBDownloadOutputs) -> None:
    for item in files.downloads:
        if item.path != files.version:
            log.info(f"Removing {item.path}")
            item.path.unlink()
