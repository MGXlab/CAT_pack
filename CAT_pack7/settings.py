"""
The internal settings classes. And subclasses in which files can be stored.
DO NOT SET END USER EDITABLE DEFAULTS HERE! only "internal defaults" are allowed (not preferred)

Validation will set these objects after validation has been done. Therefore,
settings.py will not import any validation methods or functions
"""
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path


@dataclass(frozen=True)
class ExecutionSettings:
    threads: int
    quiet: bool
    verbose: bool
    debug: bool


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
            "--no-self-hits" if self.no_self_hits else ""
        ]
        return [str(arg) for arg in command if arg != ""]


@dataclass(frozen=True)
class MMseqsSettings:
    pass


@dataclass(frozen=True)
class AnnotationFiles:
    contigs: Path
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
    contig_report: Path


@dataclass(frozen=True)
class BatFiles(AnnotationFiles):
    bin_report: Path
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


# Database preparation settings and outputs.
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
    diamond: Path | None
    cleanup: bool
