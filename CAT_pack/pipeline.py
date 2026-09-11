from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from utils.errors import ExternalToolError, InputError
from typing import Protocol
from time import sleep


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
    proteins: Path
    # some of the optional args
    range_: Decimal
    output_prefix: Path
    threads: int = 1

# Just placeholders for now. Stagebuilding will come later
# global for now to allow for easier import into cli.py
stages = ["validate_input", 'protein_prediciton',
              "alignment", 'classify']


def run_cat(aruments: CatArgs, report: Report) -> str:



    for index, stage in enumerate(stages):
        if stage == "validate_input":
            for i in range(100):
                report(stage, "running", i+1, 100)
                sleep(0.1)




        if stage == "protein_prediciton" and aruments.proteins is not None:
            report(stage, "file reused", 0 , None)

    return "done"