import logging
import shutil
from itertools import islice
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from CAT_pack.io.parsers import FastaParser, configure_polars_threads, optionally_compressed_handle
from . import tax
from .config.settings import PrepareSettings
from .tools.tool_runner import run_tool
from .utils.check import check_outputs
from .utils.errors import InputError
from .utils.logging import Report, Status

log = logging.getLogger("CAT_pack")

ACCESSION_BATCH_ROWS = 1_000_000
LCA_RECORD_BATCH = 100_000


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


def check_taxonomy_cycles(taxonomy: tax.Taxonomy, path: Path) -> None:
    # state: 
    #    0 = unvisited
    #    1 = on the current walk
    #    2 = known to reach the root
    parents = taxonomy.parents
    state = np.zeros(len(parents), dtype=np.uint8)
    for start in range(len(parents)):
        if state[start] != 0:
            continue
        node = start
        while state[node] == 0:
            state[node] = 1
            parent = int(parents[node])
            if parent == node:
                break
            node = parent
        if state[node] == 1 and int(parents[node]) != node:
            raise InputError(
                f"Cycle in taxonomy at taxid {taxonomy.taxids[node]}.",
                path=path,
                hint="Please make sure that the input files are valid.",
            )
        node = start
        while state[node] == 1:
            state[node] = 2
            parent = int(parents[node])
            if parent == node:
                break
            node = parent


def import_fasta_headers(
    fasta_file: Path, scratch: Path, report: Report, batch_rows: int = ACCESSION_BATCH_ROWS,
) -> int:
    """Scan FASTA headers into Parquet batches under scratch; return header count."""
    import polars as pl

    fastaid2prot_accessions_dir = scratch / "fastaid2prot_accessions"
    fastaid_dir = scratch / "fastaid"
    fastaid2prot_accessions_dir.mkdir()
    fastaid_dir.mkdir()

    log.info(f"Loading file {fasta_file}.")
    fastaid2prot_accessions: set[str] = set()
    records: list[int] = []
    positions: list[int] = []
    prot_accessions_batch: list[str] = []
    id_records: list[int] = []
    fastaids: list[str] = []
    total = 0
    part = 0

    def flush() -> None:
        nonlocal part, records, positions, prot_accessions_batch, id_records, fastaids
        if not prot_accessions_batch and not fastaids:
            return
        if prot_accessions_batch:
            pl.DataFrame({
                "record": records,
                "pos": positions,
                "accession": prot_accessions_batch,
            }).write_parquet(fastaid2prot_accessions_dir / f"part_{part:05d}.parquet")
        if fastaids:
            pl.DataFrame({
                "record": id_records,
                "fastaid": fastaids,
            }).write_parquet(fastaid_dir / f"part_{part:05d}.parquet")
        part += 1
        records, positions, prot_accessions_batch = [], [], []
        id_records, fastaids = [], []

    for prot_accessions in FastaParser(fasta_file).iter_header_accessions():
        fastaid = prot_accessions[0]
        if fastaid in fastaid2prot_accessions:
            raise InputError(f"Duplicate database FASTA identifier: {fastaid}.", path=fasta_file)
        fastaid2prot_accessions.add(fastaid)
        for pos, prot_accession in enumerate(prot_accessions):
            records.append(total)
            positions.append(pos)
            prot_accessions_batch.append(prot_accession)
        id_records.append(total)
        fastaids.append(fastaid)
        total += 1
        if total % 10000 == 0:
            report("Make fastaid2LCAtaxid", Status.RUNNING, total, None)
        if len(prot_accessions_batch) >= batch_rows:
            flush()

    flush()
    if total == 0:
        raise InputError("Database FASTA contains no headers.", path=fasta_file)
    return total


def get_accession_columns(prot_accession2taxid_file: Path) -> tuple[int, int]:
    with optionally_compressed_handle(prot_accession2taxid_file, "rt") as f1:
        columns = f1.readline().rstrip("\r\n").split("\t")
    try:
        return columns.index("accession.version"), columns.index("taxid")
    except ValueError as error:
        raise InputError(
            "Accession table needs accession.version and taxid columns.",
            path=prot_accession2taxid_file,
        ) from error


def bz2_to_parquet(
    prot_accession2taxid_file: Path, scratch: Path, report: Report,
) -> Path:
    """Stage a .bz2 accession table as Parquet parts, because 
    polars can't scan bz2 directly."""
    import polars as pl

    accession_column, taxid_column = get_accession_columns(prot_accession2taxid_file)
    selected = [f"field_{accession_column}", f"field_{taxid_column}"]
    max_column = max(accession_column, taxid_column)
    column_names = [f"field_{i}" for i in range(max_column + 1)]
    schema = dict.fromkeys(column_names, pl.String)
    reader_options = dict(
        has_header=False, separator="\t", quote_char=None, schema=schema,
        empty_string_is_null=False, raise_if_empty=False,
        missing_columns="insert", extra_columns="ignore",
    )
    out_dir = scratch / "prot_accession2taxid"
    out_dir.mkdir()
    part = 0
    rows = 0
    chunk_rows = 500_000

    def write_chunk(lines: list[str]) -> None:
        nonlocal part, rows
        frame = (
            pl.read_csv("".join(lines).encode("utf-8"), **reader_options)
            .select(selected)
            .rename({selected[0]: "accession.version", selected[1]: "taxid"})
            .drop_nulls(["accession.version", "taxid"])
        )
        if frame.height:
            frame.write_parquet(out_dir / f"part_{part:05d}.parquet")
            part += 1
            rows += frame.height

    with optionally_compressed_handle(prot_accession2taxid_file, "rt") as f1:
        f1.readline()
        lines = (line for line in f1 if line.strip())
        number = 0
        while chunk := list(islice(lines, chunk_rows)):
            number += 1
            write_chunk(chunk)
            report("Make fastaid2LCAtaxid", Status.RUNNING, number * chunk_rows, None)
    log.info(f"Loaded {rows:,} accession table rows from {prot_accession2taxid_file}.")
    return out_dir


def import_prot_accession2taxid(
    prot_accession2taxid_file: Path, scratch: Path, report: Report,
):
    import polars as pl

    log.info(f"Loading file {prot_accession2taxid_file}.")
    if prot_accession2taxid_file.suffix == ".bz2":
        prot_accession2taxid = pl.scan_parquet(
            str(bz2_to_parquet(prot_accession2taxid_file, scratch, report) / "*.parquet"),
        )
    else:
        get_accession_columns(prot_accession2taxid_file)  # validates the header
        prot_accession2taxid = (
            pl.scan_csv(
                prot_accession2taxid_file,
                glob=False,
                has_header=True,
                separator="\t",
                quote_char=None,
                infer_schema=False,
                truncate_ragged_lines=True,
            )
            .select(["accession.version", "taxid"])
            .drop_nulls(["accession.version", "taxid"])
        )
    return prot_accession2taxid.unique(subset=["accession.version"], keep="last")


def get_taxids(
    scratch: Path, prot_accession2taxid, *, primary_only: bool = False,
):
    import polars as pl

    fastaid2prot_accessions = pl.scan_parquet(
        str(scratch / "fastaid2prot_accessions" / "*.parquet"),
    )
    if primary_only:
        fastaid2prot_accessions = (
            fastaid2prot_accessions
            .filter(pl.col("pos") == 0)
            .select(["record", "accession"])
        )
        return (
            fastaid2prot_accessions.join(
                prot_accession2taxid,
                left_on="accession",
                right_on="accession.version",
                how="left",
            )
            .select(["record", pl.col("taxid").alias("primary_taxid")])
            .collect(engine="streaming")
        )
    return (
        fastaid2prot_accessions.join(
            prot_accession2taxid,
            left_on="accession",
            right_on="accession.version",
            how="inner",
        )
        .select(["record", "pos", "taxid"])
        .collect(engine="streaming")
    )


def taxids_to_index(joined, taxonomy: tax.Taxonomy):
    import polars as pl

    taxid_frame = pl.DataFrame({
        "taxid": taxonomy.taxids,
        "taxid_i": np.arange(len(taxonomy.taxids), dtype=np.int32),
    })
    return (
        joined.lazy()
        .join(taxid_frame.lazy(), on="taxid", how="inner")
        .select(["record", "pos", "taxid_i"])
        .sort(["record", "pos"])
        .collect(engine="streaming")
    )


def find_LCA_batch(
    taxonomy: tax.Taxonomy,
    record_ids: np.ndarray,
    taxid_indices: np.ndarray,
    record_start: int,
    record_end: int,
) -> tuple[np.ndarray, int | None]:
    batch_size = record_end - record_start
    if batch_size == 0:
        return np.empty(0, dtype=np.int32), None
    if len(record_ids) == 0:
        return np.full(batch_size, taxonomy.no_taxid, dtype=np.int32), None

    local_records = record_ids - record_start
    counts = np.bincount(local_records, minlength=batch_size)
    offsets = np.empty(batch_size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    ordered = np.argsort(local_records, kind="stable")
    hits = taxid_indices[ordered]
    LCAtaxids, _ = taxonomy.find_LCA_for_ORF(offsets, hits)
    disconnected = (counts > 0) & (LCAtaxids == taxonomy.no_taxid)
    if np.any(disconnected):
        bad = int(record_start + np.flatnonzero(disconnected)[0])
        return LCAtaxids, bad
    return LCAtaxids, None


def make_fastaid2LCAtaxid_file(
    settings: PrepareSettings, report: Report,
) -> None:
    configure_polars_threads(settings.threads)
    import polars as pl

    output = settings.files.fastaid2LCAtaxid
    scratch_parent = output.parent
    scratch_parent.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(prefix=".CAT_pack_fastaid2LCA_", dir=scratch_parent) as scratch_name:
        scratch = Path(scratch_name)
        total = import_fasta_headers(settings.db_fasta, scratch, report, ACCESSION_BATCH_ROWS)
        taxonomy = tax.Taxonomy.import_nodes(settings.files.nodes)
        check_taxonomy_cycles(taxonomy, settings.files.nodes)

        log.info("Finding LCA of all protein accession numbers in fasta headers.")
        prot_accession2taxid_table = import_prot_accession2taxid(settings.acc2tax, scratch, report)
        joined = get_taxids(scratch, prot_accession2taxid_table)
        primary_taxids = get_taxids(
            scratch, prot_accession2taxid_table, primary_only=True,
        )
        hits = taxids_to_index(joined, taxonomy)
        del joined

        hit_records = hits.get_column("record").to_numpy()
        hit_taxids = hits.get_column("taxid_i").to_numpy()
        del hits

        fastaid2record = (
            pl.scan_parquet(str(scratch / "fastaid" / "*.parquet"))
            .sort("record")
            .collect(engine="streaming")
        )
        fastaids = fastaid2record.get_column("fastaid").to_list()
        del fastaid2record

        # Taxid of the primary (first) accession of each header, keyed by fastaid.
        prot_accession2taxid: dict[str, str] = {}
        for record, taxid in zip(
            primary_taxids.get_column("record").to_list(),
            primary_taxids.get_column("primary_taxid").to_list(),
        ):
            if taxid is not None:
                prot_accession2taxid[fastaids[int(record)]] = taxid
        del primary_taxids

        no_taxid = 0
        corrected = 0
        with output.open("w", encoding="utf-8", newline="\n") as outf1:
            for batch_start in range(0, total, LCA_RECORD_BATCH):
                batch_end = min(batch_start + LCA_RECORD_BATCH, total)
                # hit_records is sorted, so a batch is a contiguous slice
                start, end = np.searchsorted(hit_records, [batch_start, batch_end])
                batch_record_ids = hit_records[start:end]
                batch_taxids = hit_taxids[start:end]
                LCAtaxids, bad_record = find_LCA_batch(
                    taxonomy, batch_record_ids, batch_taxids, batch_start, batch_end,
                )
                if bad_record is not None:
                    raise InputError(
                        f"Cannot find a shared taxonomic ancestor for the protein "
                        f"accessions keys in FASTA header {fastaids[bad_record]}.",
                        path=settings.files.nodes,
                        hint="Their taxonomic lineages should share a root taxid. "
                             "please use a valid input file",
                    )
                for offset, lca_index in enumerate(LCAtaxids):
                    record = batch_start + offset
                    if (record + 1) % 10000 == 0 or record + 1 == total:
                        report("Make fastaid2LCAtaxid", Status.RUNNING, record + 1, total)
                    if lca_index == taxonomy.no_taxid:
                        no_taxid += 1
                        continue
                    LCAtaxid = taxonomy.taxids[int(lca_index)]
                    fastaid = fastaids[record]
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
    settings: PrepareSettings, report: Report,
) -> None:
    taxid2parent = tax.read_nodes(settings.files.nodes)
    taxid2offspring = find_offspring(settings.files.fastaid2LCAtaxid, taxid2parent, report)
    taxids_with_multiple_offspring_file = settings.files.taxids_with_multiple_offspring
    log.info(f"Writing {taxids_with_multiple_offspring_file}.")
    with taxids_with_multiple_offspring_file.open("w", encoding="utf-8", newline="\n") as outf1:
        for taxid in taxid2offspring:
            if len(taxid2offspring[taxid]) >= 2:
                outf1.write(f"{taxid}\n")


def make_diamond_database(settings: PrepareSettings, report: Report) -> None:
    log.info(f"Constructing DIAMOND database {settings.files.diamond_database} "
             f"from {settings.db_fasta} using {settings.threads} cores.")
    run_tool(settings.diamond.get_command(), "DIAMOND")


def make_mmseqs2_database(settings: PrepareSettings, report: Report) -> None:
    log.info(f"Constructing MMseqs2 database {settings.files.mmseqs2_database} "
             f"from {settings.db_fasta} using {settings.threads} cores.")
    run_tool(settings.mmseqs.get_command(), "MMseqs2")


def make_mmseqs2_gpu_database(settings: PrepareSettings, report: Report) -> None:
    database = settings.files.mmseqs2_gpu_database.resolve()
    log.info("Constructing GPU MMseqs2 database %s from %s.",
             database, settings.files.mmseqs2_database)
    with TemporaryDirectory(prefix=".CAT_pack_MMseqs2_GPU_", dir=database.parent) as scratch:
        temporary_database = Path(scratch) / "db"
        run_tool(settings.mmseqs.get_gpu_command(temporary_database), "MMseqs2",
                 operation="GPU database preparation")
        outputs = sorted(Path(scratch).iterdir(), key=lambda path: path.name == "db.dbtype")
        destinations = [database.with_name(database.name + path.name[2:]) for path in outputs]
        check_outputs(destinations)
        for source, destination in zip(outputs, destinations):
            source.replace(destination)


def copy_taxonomy(settings: PrepareSettings) -> None:
    for source, target in ((settings.names, settings.files.names), (settings.nodes, settings.files.nodes)):
        if target.is_file():
            continue
        shutil.copyfile(source, target)


