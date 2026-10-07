"""
Here will al the user configurable (with default) classes live, default values can be set here,
and within cli.py classname.value can be called to set the default value correctly
for the argument parser. So that you only have to change the defaults here.
And can easily build or extend classes.
"""
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path
from typing import Literal

AlignerName = Literal["diamond", "mmseqs2"]


@dataclass(frozen=True, kw_only=True)
class ExecutionOptions:
    quiet: bool = False
    verbose: bool = False
    debug: bool = False
    threads: int = 1


@dataclass(frozen=True)
class DiamondOptions:
    mode: str = "default"
    no_self_hits: bool = False
    block_size: float = 12.0
    index_chunks: int = 1
    path_to_diamond: Path | None = None


@dataclass(frozen=True)
class MMseqsOptions:
    sensitivity: float = 5.7
    split_memory_limit: str = "0"
    executable: Path | None = None


@dataclass(frozen=True)
class MemoryOptions:
    low_memory: bool = False
    high_memory: bool = False
    available_memory: int | None = None


@dataclass(frozen=True, kw_only=True)
class AnnotationOptions(ExecutionOptions):
    database: Path
    proteins: Path | None = None
    alignment: Path | None = None
    range_: Decimal = Decimal("10.0")
    fraction: Decimal = Decimal("0.5")
    log_file: Path | None = None
    output_prefix: Path = Path("out.CAT")
    diamond: DiamondOptions = field(default_factory=DiamondOptions)
    mmseqs: MMseqsOptions = field(default_factory=MMseqsOptions)
    aligner: AlignerName = "diamond"
    top: int = 11
    tmpdir: Path | None = None
    compress: bool = False
    memory: MemoryOptions = field(default_factory=MemoryOptions)

    @property
    def log_path(self) -> Path:
        return self.log_file if self.log_file is not None else Path(f"{self.output_prefix}.log")


@dataclass(frozen=True, kw_only=True)
class CatOptions(AnnotationOptions):
    contigs: Path


@dataclass(frozen=True, kw_only=True)
class BatOptions(AnnotationOptions):
    bins: Path
    bin_suffix: str = ".fna"
    range_: Decimal = Decimal("5")
    fraction: Decimal = Decimal("0.3")
    output_prefix: Path = Path("out.BAT")
    no_stars: bool = False


@dataclass(frozen=True, kw_only=True)
class PrepareOptions(ExecutionOptions):
    db_fasta: Path
    names: Path
    nodes: Path
    acc2tax: Path
    db_dir: Path
    path_to_diamond: Path | None = DiamondOptions.path_to_diamond
    common_prefix: str | None = None
    build_mmseqs2: bool = False
    path_to_mmseqs: Path | None = None
