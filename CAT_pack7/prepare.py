import gzip
import logging
import shutil
from pathlib import Path
from typing import TextIO

from . import tax
from .settings import PrepareSettings
from .tools.tool_runner import run_tool
from .utils.errors import InputError
from .utils.logging import Report, Status

log = logging.getLogger("CAT_pack")


def optionally_compressed_handle(file_path: Path) -> TextIO:
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt", encoding="utf-8")
    return file_path.open(encoding="utf-8")


def find_lineage(taxid: str, taxid2parent: dict[str, str]) -> list[str]:
    lineage: list[str] = []
    flown_by_taxids: set[str] = set()
    while True:
        if taxid in flown_by_taxids:
            raise InputError(f"Cycle in taxonomy at taxid {taxid}.",
                             hint="Please make sure that the input files are valid.")
        flown_by_taxids.add(taxid)
        lineage.append(taxid)
        parent = taxid2parent[taxid]
        if parent == taxid:
            return lineage
        taxid = parent


def import_fasta_headers(fasta_file: Path, report: Report) -> tuple[dict[str, list[str]], set[str]]:
    log.info(f"Loading file {fasta_file}.")
    fastaid2prot_accessions: dict[str, list[str]] = {}
    prot_accessions_whitelist: set[str] = set()
    with optionally_compressed_handle(fasta_file) as f1:
        for line in f1:
            if not line.startswith(">"):
                continue
            prot_accessions: list[str] = []
            for part in line[1:].strip().split("\x01"):
                fields = part.split()
                if not fields:
                    raise InputError("Database FASTA contains an empty header.", path=fasta_file)
                prot_accessions.append(fields[0])
            fastaid = prot_accessions[0]
            if fastaid in fastaid2prot_accessions:
                raise InputError(f"Duplicate database FASTA identifier: {fastaid}.", path=fasta_file)
            fastaid2prot_accessions[fastaid] = prot_accessions
            prot_accessions_whitelist.update(prot_accessions)
            if len(fastaid2prot_accessions) % 10000 == 0:
                report("Make fastaid2LCAtaxid", Status.RUNNING, len(fastaid2prot_accessions), None)
    if not fastaid2prot_accessions:
        raise InputError("Database FASTA contains no headers.", path=fasta_file)
    return fastaid2prot_accessions, prot_accessions_whitelist


def import_prot_accession2taxid(
    prot_accession2taxid_file: Path, prot_accessions_whitelist: set[str], report: Report
) -> dict[str, str]:
    log.info(f"Loading file {prot_accession2taxid_file}.")
    prot_accession2taxid: dict[str, str] = {}
    with optionally_compressed_handle(prot_accession2taxid_file) as f1:
        columns = f1.readline().rstrip().split("\t")
        try:
            accession_column = columns.index("accession.version")
            taxid_column = columns.index("taxid")
        except ValueError as error:
            raise InputError("Accession table needs accession.version and taxid columns.", path=prot_accession2taxid_file) from error
        for number, line in enumerate(f1, 2):
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) <= max(accession_column, taxid_column):
                raise InputError(f"Incomplete accession table row {number}.", path=prot_accession2taxid_file)
            prot_accession = fields[accession_column]
            if prot_accession in prot_accessions_whitelist:
                prot_accession2taxid[prot_accession] = fields[taxid_column]
            if number % 100000 == 0:
                report("Make fastaid2LCAtaxid", Status.RUNNING, number, None)
    return prot_accession2taxid


def make_fastaid2LCAtaxid_file(
    fastaid2LCAtaxid_file: Path, fasta_file: Path, prot_accession2taxid_file: Path,
    taxid2parent: dict[str, str], report: Report
) -> None:
    fastaid2prot_accessions, prot_accessions_whitelist = import_fasta_headers(fasta_file, report)
    prot_accession2taxid = import_prot_accession2taxid(
        prot_accession2taxid_file, prot_accessions_whitelist, report,
    )
    log.info("Finding LCA of all protein accession numbers in fasta headers.")
    no_taxid = 0
    corrected = 0
    total = 0
    with fastaid2LCAtaxid_file.open("w", encoding="utf-8", newline="\n") as outf1:
        for fastaid, prot_accessions in fastaid2prot_accessions.items():
            list_of_lineages: list[list[str]] = []
            for prot_accession in prot_accessions:
                try:
                    taxid = prot_accession2taxid[prot_accession]
                    lineage = find_lineage(taxid, taxid2parent)
                    list_of_lineages.append(lineage)
                except KeyError:
                    continue
            total += 1
            if total % 10000 == 0:
                report("Make fastaid2LCAtaxid", Status.RUNNING, total, len(fastaid2prot_accessions))
            if len(list_of_lineages) == 0:
                no_taxid += 1
                continue

            LCAtaxid = tax.find_LCA(list_of_lineages)
            if LCAtaxid is None:
                raise InputError(f"Disconnected taxonomic lineages in FASTA header {fastaid}.", path=fasta_file)
            outf1.write(f"{fastaid}\t{LCAtaxid}\n")
            if fastaid not in prot_accession2taxid or LCAtaxid != prot_accession2taxid[fastaid]:
                corrected += 1
    log.info(f"Mapped {total - no_taxid:,}/{total:,} headers; "
             f"{corrected:,} corrected using secondary accessions; {no_taxid:,} without usable taxonomy.")


def find_offspring(fastaid2LCAtaxid_file: Path, taxid2parent: dict[str, str], report: Report) -> dict[str, set[str]]:
    log.info("Searching database for taxids with multiple offspring.")
    taxid2offspring: dict[str, set[str]] = {}
    visited_taxids: set[str] = set()
    with fastaid2LCAtaxid_file.open(encoding="utf-8") as f1:
        for number, line in enumerate(f1, 1):
            fields = line.rstrip().split("\t")
            if len(fields) != 2:
                raise InputError(f"Invalid protein-to-taxid mapping row {number}.", path=fastaid2LCAtaxid_file)
            taxid = fields[1]
            if taxid in visited_taxids:
                continue
            visited_taxids.add(taxid)
            try:
                lineage = find_lineage(taxid, taxid2parent)
            except KeyError as error:
                raise InputError(f"Mapping taxid {taxid} has missing ancestry.", path=fastaid2LCAtaxid_file) from error
            for i, taxid in enumerate(lineage):
                if i == 0:
                    continue
                if taxid not in taxid2offspring:
                    taxid2offspring[taxid] = set()
                offspring = lineage[i - 1]
                taxid2offspring[taxid].add(offspring)
            if number % 10000 == 0:
                report("Make taxids with multiple offspring", Status.RUNNING, number, None)
    return taxid2offspring


def write_taxids_with_multiple_offspring_file(
    taxids_with_multiple_offspring_file: Path, taxid2offspring: dict[str, set[str]],
) -> None:
    log.info(f"Writing {taxids_with_multiple_offspring_file}.")
    with taxids_with_multiple_offspring_file.open("w", encoding="utf-8", newline="\n") as outf1:
        for taxid in taxid2offspring:
            if len(taxid2offspring[taxid]) >= 2:
                outf1.write(f"{taxid}\n")


def make_diamond_database(settings: PrepareSettings) -> None:
    log.info(f"Constructing DIAMOND database {settings.files.diamond_database} "
             f"from {settings.db_fasta} using {settings.threads} cores.")
    run_tool(settings.diamond.get_command(), "DIAMOND")


def make_mmseqs2_database(settings: PrepareSettings) -> None:
    log.info(f"Constructing MMseqs2 database {settings.files.mmseqs2_database} "
             f"from {settings.db_fasta} using {settings.threads} cores.")
    run_tool(settings.mmseqs.get_command(), "MMseqs2")


def copy_taxonomy(settings: PrepareSettings) -> None:
    for source, target in ((settings.names, settings.files.names), (settings.nodes, settings.files.nodes)):
        if target.is_file():
            continue
        shutil.copyfile(source, target)


