import os
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Protocol
from time import sleep

from utils.errors import ValidationError
from validation import validate_args


class Report(Protocol):
    """Reports back to cli.py with the current progress.
    Call it with a step name, the current status, how much is completed,
    how much total work there must be done (including completed) Leave None if
    the amount of work is not known (yet).

    Within status make the choice: Running, complete, reused, skipped or failed
    Cancelled can also be used if user canceld the run
    """

    def __call__(self, stage: str, status: str, completed: int,
                 total: int | None) -> None: ...


@dataclass(frozen=True)
class CatArgs:
    """Arguments for the CAT run, provided (for now) only by cli.py
    Does not validate, that happens within run_cat right know
    @bastiaan, maybe the validation can move here in the future?
    But might become messy
    """

    contigs: Path
    database: Path
    taxonomy: Path
    proteins: Path | None
    alignment: Path | None
    range_: Decimal
    fraction: Decimal
    log_file: Path
    output_prefix: Path
    threads: int = 1

@dataclass(frozen=True)
class Step:
    """
    A simple dataclass currently acting like a dictionairy that stores the
    name of the step and reuse, if reuse is true it will use userprovided data
    """

    name: str
    reuse: bool = False

    def __str__(self) -> str:
        return self.name


def build_plan(args: CatArgs) -> list[Step]:
    return [
        Step("Input validation"),
        Step("Protein prediction", reuse=args.proteins is not None),
        Step("Alignment", reuse=args.alignment is not None),
        Step("Classify"),
    ]

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
