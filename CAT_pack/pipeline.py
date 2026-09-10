from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CatArgs:
    #all the required
    contigs: Path
    database: Path
    taxonomy: Path
    # some of the optional args
    output_prefix: Path
    threads: int = 1