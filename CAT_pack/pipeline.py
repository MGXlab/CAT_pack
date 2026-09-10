from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path


@dataclass(frozen=True)
class CatArgs:
    #all the required
    contigs: Path
    database: Path
    taxonomy: Path
    # some of the optional args
    _range: Decimal
    output_prefix: Path
    threads: int = 1