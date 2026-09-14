import os
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Protocol
import traceback

from .utils.errors import CatError, InputError
from .validation import validate_args
from .tools.pyrodigal import run_pyrodigal
from .tools.diamond import run_diamond
from .classification import contig_classification


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


def run_cat(args: CatArgs, report: Report) -> dict[str, Path]:
    """Contig annotation tool (CAT) run"""

    plan = build_plan(args)
    step_index = 0
    current_step = plan[step_index]
    log = None #only activate log after input validation :D


    try:
        report(current_step.name, "running", 0, 1)
        files = validate_args(args)

        # Exclusive creation ("x") this ensures no old log overwrite

        log = files["log"].open("x", encoding="utf-8")
        log.write(f"CAT_pack7\n{args!r}\n")
        report(current_step.name, "complete", 1, 1)

        for step_index, current_step in enumerate(plan[1:], start=1):
            if current_step.reuse:
                # This forces "completion" of reused progressbar.
                # And thus makes the progressbar green, currently if only
                # reused is called the bar says dimmed (also a good indicator)
                report(current_step.name, "running", 1, 1)
                report(current_step.name, "reused", 1, 1)
                log.write(f"Reused: {current_step.name}\n")
                log.flush()
                continue

            report(current_step.name, "running", 0, None)
            log.write(f"Starting: {current_step.name}\n")
            log.flush()

            if current_step.name == "Protein prediction":
                run_pyrodigal(args, files, log, report)
            elif current_step.name == "Alignment":
                run_diamond(args, files, log, report)
            elif current_step.name == "Classify":
                contig_classification(args, files, log, report)
            else:
                log.write(f"Can't spell that well, my misspelling:"
                          f"{current_step.name}\n")

            # always complete the step with 1 of 1 so 100% completion,
            # errors will have caused the code to stop before this if there is
            # an exception.
            report(current_step.name, "complete", 1, 1)
            log.write(f"Completed: {current_step.name}\n")
            log.flush()

        return {
            "Contig classifications": Path("Future path to C2C file :)"),
            "ORF classifications": Path("Future path to ORF2LCA file :)"),
            "Log": files["log"],
        }

    except KeyboardInterrupt:
        report(current_step.name, "cancelled", 0, None)

        for skipped_steps in plan[step_index + 1:]:
            report(skipped_steps.name, "skipped", 0, None)

        if log is not None: log.write("Run cancelled\n")
        raise

    except Exception as error:
        report(current_step.name, "failed", 0, None)
        for remaining_step in plan[step_index + 1:]:
            report(remaining_step.name, "skipped", 0, None)

        # Log all details into the log file, and the CLI shows the short error
        if log is not None: traceback.print_exc(file=log)

        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = args.log_file if log is not None else None
            raise

        # Catch all other exceptions
        raise

    finally:
        if log is not None:
            log.close()
