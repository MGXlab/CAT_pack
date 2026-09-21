"""
Here will al the "Default" classes live, default values can be set here,
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
class Defaults:
    quiet: bool = False
    verbose: bool = False
    debug: bool = False
    threads: int = 1


@dataclass(frozen=True)
class DiamondDefaults:
    mode: str = "default"
    no_self_hits: bool = False
    block_size: float = 12.0
    index_chunks: int = 1
    path_to_diamond: Path | None = None


@dataclass(frozen=True)
class MMseqsDefaults:
    sensitivity: float = 5.7
    split_memory_limit: str = "0"
    executable: Path | None = None


@dataclass(frozen=True, kw_only=True)
class CatDefaults(Defaults):
    contigs: Path
    database: Path
    proteins: Path | None = None
    alignment: Path | None = None
    range_: Decimal = Decimal("10.0")
    fraction: Decimal = Decimal("0.5")
    log_file: Path | None = None
    output_prefix: Path = Path("out.CAT")
    diamond: DiamondDefaults = field(default=DiamondDefaults)
    mmseqs: MMseqsDefaults = field(default=MMseqsDefaults)
    aligner: AlignerName = "diamond"
    top: int = 11
    tmpdir: Path | None = None
    compress: bool = False

    @property
    def log_path(self) -> Path:
        return self.log_file if self.log_file is not None else Path(f"{self.output_prefix}.log")


@dataclass(frozen=True, kw_only=True)
class PrepareDefaults(Defaults):
    db_fasta: Path
    names: Path
    nodes: Path
    acc2tax: Path
    db_dir: Path
    path_to_diamond: Path | None = DiamondDefaults.path_to_diamond
    common_prefix: str | None = None
    cleanup: bool = False
