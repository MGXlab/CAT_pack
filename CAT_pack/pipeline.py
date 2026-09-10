from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from utils.errors import ExternalToolError, InputError
from typing import Protocol


class Report(Protocol):
    def __call__(
        self,
        stage: str,
        status: str,
        completed: int,
        total: int | None,
    ) -> None: ...


@dataclass(frozen=True)
class CatArgs:
    #all the required
    contigs: Path
    database: Path
    taxonomy: Path
    # some of the optional args
    range_: Decimal
    output_prefix: Path
    threads: int = 1




def run_cat(aruments: CatArgs, report: Report) -> None:

    # Just placeholders for now. Stagebuilding will come later
    stages = ["validate_input", 'protein_prediciton',
              "alignment", 'classify']

    for index, stage in enumerate(stages):
        if stage == "protein_prediciton" and aruments.proteins is not None:
            report(stage, "file reused", 0 , None)