import logging
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from .classification import contig_classification
from .tools.aligner import run_aligner, DiamondArgs, MMseqsArgs
from .tools.pyrodigal import run_protein_prediction
from .utils.errors import CatError
from .utils.logging import Status, Report
from .validation import validate_args, validate_aligner_args


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
    aligner: DiamondArgs | MMseqsArgs | None = None
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
    log = logging.getLogger("CAT_pack")


    try:
        report(current_step.name, Status.RUNNING, 0, 1)
        files = validate_args(args)
        aligner_args = validate_aligner_args(args)
        # Exclusive creation ("x") this ensures no old log overwrite


        report(current_step.name, Status.COMPLETE, 1, 1)

        for step_index, current_step in enumerate(plan[1:], start=1):
            if current_step.reuse:
                # This forces "completion" of reused progressbar.
                # And thus makes the progressbar green, currently if only
                # reused is called the bar says dimmed (also a good indicator)
                report(current_step.name, Status.SUPPLIED, 1, 1)
                log.info(f"Supplied file: {current_step.name}")
                continue

            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")

            if current_step.name == "Protein prediction":
                run_protein_prediction(args, files, report, "pyrodigal")
            elif current_step.name == "Alignment":
                run_aligner(aligner_args, report)
            elif current_step.name == "Classify":
                contig_classification(args, files, report)
            else:
                log.error(f"Can't spell that well, my misspelling:"
                          f"{current_step.name}")

            # always complete the step with 1 of 1 so 100% completion,
            # errors will have caused the code to stop before this if there is
            # an exception.
            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")

        return {
            "Contig classifications": Path("Future path to C2C file :)"),
            "ORF classifications": Path("Future path to ORF2LCA file :)"),
            "Log": args.log_file if args.log_file is not None else None,
        }

    except KeyboardInterrupt:
        report(current_step.name, Status.CANCELLED, 0, None)

        for skipped_steps in plan[step_index + 1:]:
            report(skipped_steps.name, Status.SKIPPED, 0, None)

        if log is not None: log.error("Run cancelled")
        raise

    except Exception as error:
        report(current_step.name, Status.FAILED, 0, None)
        for remaining_step in plan[step_index + 1:]:
            report(remaining_step.name, Status.SKIPPED, 0, None)

        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = args.log_file if log is not None else None
            raise

        # Catch all other exceptions
        raise

