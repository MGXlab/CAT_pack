import os
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


def check_folder(path: Path, label: str) -> None:
    if not path.is_dir():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the folder exists.",
            path=path,
        )


def validate_args(args):
    check_folder(args.database, "Database folder")




def run_cat(args: CatArgs, report: Report) -> str:
    step = steps[0]

    report(step, "running", 0, 1)

    try:
        report(step, "running", 0, 1)
        validate_args(args)
        report(step, "complete", 1, 1)
    except Exception as e:
        report(step, "failed", 0, None)

        for remaining in steps[1:]:
            report(remaining, "skipped", 0, None)

        raise
    report(step, "complete", 1, 1)
