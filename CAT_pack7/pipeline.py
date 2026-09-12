import os
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Protocol
from time import sleep

from utils.errors import ValidationError
from validation import validate_args


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
    proteins: Path | None
    # some of the optional args
    range_: Decimal
    fraction: Decimal
    log_file: Path
    output_prefix: Path
    threads: int = 1

# Just placeholders for now. Stagebuilding will come later
# global for now to allow for easier import into cli.py
steps = ["Input validation", 'protein_prediciton',
              "alignment", 'classify']


def run_cat(args: CatArgs, report: Report) -> str:
    step = steps[0]

    report(step, "running", 0, 1)

    try:
        validate_args(args)
    except ValidationError:
        report(step, "failed", 0, None)

        for remaining in steps[1:]:
            report(remaining, "skipped", 0, None)

        raise
    report(step, "complete", 1, 1)
