"""
The internal settings classes. And subclasses in which files can be stored.
DO NOT SET END USER EDITABLE DEFAULTS HERE! only "internal defaults" are allowed (not preferred)

Validation will set these objects after validation has been done. Therefore,
settings.py will not import any validation methods or functions
"""
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path


@dataclass(frozen=True)
class MemorySettings:
    low_memory: bool = False
    available_memory_bytes: int | None = None


@dataclass(frozen=True)
class MemoryBudget:
    total_bytes: int
    workers: int
    batch_bytes: int


@dataclass(frozen=True)
class ExecutionSettings:
    threads: int
    quiet: bool
    verbose: bool
    debug: bool
    memory: MemorySettings = field(default_factory=MemorySettings, kw_only=True)


@dataclass(frozen=True)
class ClassificationSettings:
    range_: Decimal
    fraction: Decimal


@dataclass(frozen=True)
class TaxonomyFiles:
    names: Path
    nodes: Path


@dataclass(frozen=True)
class DatabaseFiles(TaxonomyFiles):
    fastaid2LCA: Path
    branches: Path
    diamond: Path | None


@dataclass(frozen=True)
class DiamondParameters:
    diamond: Path
    mode: str
    no_self_hits: bool
    block_size: float
    index_chunks: int


@dataclass(frozen=True)
class DiamondSettings(DiamondParameters):
    query: Path
    database: Path
    alignment: Path
    tmpdir: Path
    threads: int
    top: int
    compression: bool
    verbose: bool
    blast_flavour: str = "blastp"
    matrix: str = "BLOSUM62"
    evalue: str = "0.001"
    cat_only: bool = True

    def get_command(self):
        command = [
            self.diamond, self.blast_flavour,
            "-q", self.query,
            "-d", self.database,
            "--top", self.top,
            "--matrix", self.matrix,
            "--evalue", self.evalue,
            "-o", self.alignment,
            "-p", self.threads,
            "--block-size", self.block_size,
            "--index-chunks", self.index_chunks,
            "--tmpdir", self.tmpdir,
            "--compress", int(self.compression),
            f"--{self.mode}" if self.mode != 'default' else "",
            "--quiet" if not self.verbose else "",
            "--no-self-hits" if self.no_self_hits else "",
        ]
        if self.cat_only:
            command.extend(["--outfmt", "6", "qseqid", "sseqid", "bitscore"])
        return [str(arg) for arg in command if arg != ""]


@dataclass(frozen=True)
class MMseqsSettings:
    pass


@dataclass(frozen=True)
class AnnotationFiles:
    report: Path
    proteins_fasta: Path
    proteins_gff: Path | None
    alignment: Path
    fastaid2LCA: Path
    branches: Path
    names: Path
    nodes: Path
    orf_report: Path
    diamond_database: Path | None


@dataclass(frozen=True)
class CatFiles(AnnotationFiles):
    contigs: Path


@dataclass(frozen=True)
class BatFiles(AnnotationFiles):
    bin2contigs: dict[str, list[str]]
    bin_paths: tuple[Path, ...]


@dataclass(frozen=True)
class CatSettings(ExecutionSettings, ClassificationSettings):
    files: CatFiles
    aligner: DiamondSettings | None
    log_file: Path


@dataclass(frozen=True)
class BatSettings(ExecutionSettings, ClassificationSettings):
    files: BatFiles
    aligner: DiamondSettings | None
    log_file: Path
    no_stars: bool


@dataclass(frozen=True)
class RatFiles:
    contigs_fasta: Path
    read_files: tuple[Path, ...]
    bam_files: tuple[Path, ...]
    bin_paths: tuple[Path, ...]
    bin2contigs: dict[str, list[str]]
    contigs2classification: Path | None
    bin2classification: Path | None
    unmapped2classification: Path | None
    unclassified_fasta: Path | None
    alignment_unmapped: Path | None
    contig_abundance: Path
    complete_abundance: Path
    bin_abundance: Path | None
    read_classification: Path | None
    mapping_sam: Path | None = None
    mapping_bam: Path | None = None
    outputs: tuple[Path, ...] = ()


@dataclass(frozen=True)
class RatSettings(ExecutionSettings):
    files: RatFiles
    log_file: Path
    bwa: Path
    samtools: Path
    sort_memory_mib: int | None = None
    cat: CatSettings | None = None
    bat: BatSettings | None = None


# Database preparation settings and outputs.
@dataclass(frozen=True)
class DiamondPrepareSettings:
    diamond: Path
    db_fasta: Path
    database: Path
    threads: int
    verbose: bool

    def get_command(self) -> list[str]:
        command = [
            str(self.diamond), "makedb",
            "--in", str(self.db_fasta.resolve()),
            "-d", str(self.database.with_suffix("").resolve()),
            "-p", str(self.threads),
        ]
        if not self.verbose:
            command.append("--quiet")
        return [str(arg) for arg in command if arg != ""]


@dataclass(frozen=True)
class MMseqsPrepareSettings:
    mmseqs: Path
    db_fasta: Path
    database: Path
    threads: int
    verbose: bool

    def get_command(self) -> list[str]:
        command = [
            str(self.mmseqs), "createdb",
            str(self.db_fasta.resolve()), str(self.database.resolve()),
            "--threads", str(self.threads),
            "--compressed", "1",
        ]
        if not self.verbose:
            command.extend(["-v", "0"])
        return command


@dataclass(frozen=True)
class PrepareOutputs:
    prefix: str
    data_folder: Path
    names: Path
    nodes: Path
    log_file: Path
    diamond_database: Path
    mmseqs2_database: Path
    fastaid2LCAtaxid: Path
    taxids_with_multiple_offspring: Path


@dataclass(frozen=True)
class PrepareSettings(ExecutionSettings):
    files: PrepareOutputs
    db_fasta: Path
    names: Path
    nodes: Path
    acc2tax: Path
    diamond: DiamondPrepareSettings | None
    mmseqs: MMseqsPrepareSettings | None = None

@dataclass()
class ParserSettings:
    # Feedback from @Bastiaan
    # DEVnote to make the multiplier dependent on the significance of the bitscore.
    bit_score_multiplier = 10



@dataclass(frozen=True)
class DownloadFile:
    url: str
    path: Path
    checksum: Path | None = None


@dataclass()
class NRUrls:
    urls: list[str]

    def __init__(self):
        # HTTPS avoids waiting for FTP's final transfer acknowledgement.
        base_tax_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"
        base_blast_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA"

        tax_file_list = ["taxdump.tar.gz", "taxdump.tar.gz.md5",
                         "accession2taxid/prot.accession2taxid.FULL.gz",
                         "accession2taxid/prot.accession2taxid.FULL.gz.md5"]
        blast_file_list = ["nr.gz", "nr.gz.md5"]
        self.urls = [f"{url}/{file}" for url, files in ((base_tax_url, tax_file_list), (base_blast_url, blast_file_list)) for file in files]






@dataclass
class GTDBUrls:
    urls: list[str]

    def __init__(self):
        # Using "latest" as an entry point.
        # This needs to be checked for future versions.
        # Use the European mirror.
        parent_url = "https://data.gtdb.aau.ecogenomic.org/releases/latest"

        file_list = [
            "VERSION.txt",
            "ar53_taxonomy.tsv.gz",
            "bac120_taxonomy.tsv.gz",
            "MD5SUM.txt",
            "bac120.tree",
            "ar53.tree",
            "genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"
        ]
        self.urls = [f"{parent_url}/{file_}" for file_ in file_list]


@dataclass(frozen=True)
class DownloadOutputs(TaxonomyFiles):
    prefix: str
    data_folder: Path
    log_file: Path
    db_fasta: Path
    acc2tax: Path
    taxonomy_archive: Path | None
    downloads: tuple[DownloadFile, ...]


@dataclass(frozen=True)
class GTDBDownloadOutputs(DownloadOutputs):
    protein_archive: Path
    taxonomy_tables: tuple[Path, Path]
    trees: tuple[Path, Path]
    tree: Path
    version: Path
    checksums: Path


@dataclass(frozen=True)
class DownloadSettings:
    db_type: str
    files: DownloadOutputs
    cleanup: bool
    quiet: bool
    no_log: bool
    debug: bool
