import logging
from dataclasses import dataclass
from pathlib import Path

from .classification import contig_classification
from .defaults import CatDefaults, PrepareDefaults
from .tools.aligner import run_aligner
from .tools.pyrodigal import run_protein_prediction
from .utils.errors import CatError
from .utils.logging import Status, Report
from .validation import get_file_names, get_validated_settings, validate_prepare


@dataclass()
class Step:
    """
    A simple dataclass currently acting like a dictionairy that stores the
    name of the step and reuse, if reuse is true it will use userprovided data
    """

    name: str
    supplied: bool = False
    citation_message: bool = False


    def __str__(self) -> str:
        return self.name


def build_plan(args: CatDefaults | PrepareDefaults, files=None) -> list[Step]:
    if type(args) == PrepareDefaults:
        return [
            Step("Input validation"),
            Step("Make DIAMOND database",
                 supplied=files is not None and files.diamond_database.is_file()),
            Step("Make MMseqs2 database",
                 supplied=files is not None and files.mmseqs2_database.is_file()),
            Step("Make fastaid2LCAtaxid",
                 supplied=files is not None and files.fastaid2LCAtaxid.is_file()),
            Step("Make taxids with multiple offspring",
                 supplied=files is not None and files.taxids_with_multiple_offspring.is_file()),
        ]
    if type(args) == CatDefaults:
        return [
            Step("Input validation"),
            Step("Protein prediction", supplied=args.proteins is not None),
            Step("Alignment", supplied=args.alignment is not None),
            Step("Classify"),
        ]
    raise CatError("I haven't figured out how to build that specific plan")


def run_cat(args: CatDefaults, report: Report) -> dict[str, Path]:
    """Contig annotation tool (CAT) run"""

    plan = build_plan(args)
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")


    try:
        report(current_step.name, Status.RUNNING, 0, 1)
        log.info(f"Starting: {current_step.name}")
        settings = get_validated_settings(args)
        report(current_step.name, Status.COMPLETE, 1, 1)
        log.info(f"Completed: {current_step.name}")

        for step_index, current_step in enumerate(plan[1:], start=1):
            if current_step.supplied:
                report(current_step.name, Status.SUPPLIED, 1, 1)
                log.info(f"Supplied file: {current_step.name}")
                continue

            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")

            if current_step.name == "Protein prediction":
                run_protein_prediction(settings, report, "pyrodigal")
            elif current_step.name == "Alignment":
                run_aligner(settings.aligner, report)
            elif current_step.name == "Classify":
                contig_classification(settings, settings.files, report)
            else:
                log.error(f"Can't spell that well, my misspelling:"
                          f"{current_step.name}")

            # always complete the step with 1 of 1 so 100% completion,
            # errors will have caused the code to stop before this if there is
            # an exception.
            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")

        return {
            "Contig classifications": settings.files.contig_report,
            "ORF classifications": settings.files.orf_report,
            "Log": settings.log_file,
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
            error.log_file = args.log_path if args.log_path.is_file() else None
            raise

        # Catch all other exceptions
        raise





def run_prepare(args: PrepareDefaults, report: Report):
    files = get_file_names(args)
    plan = build_plan(args, get_file_names(args))
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")

    try:
        for step_index, current_step in enumerate(plan):
            if current_step.supplied:
                report(current_step.name, Status.SUPPLIED, 1, 1)
                log.info(f"Already exists, skipped making of: {current_step.name}")
                continue

            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")
            if step_index == 0:
                settings = validate_prepare(args)
            # TODO: Port over the next steps
            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")

        return 0 # so succes :)

    except KeyboardInterrupt:
        report(current_step.name, Status.CANCELLED, 0, None)
        for skipped_steps in plan[step_index + 1:]:
            report(skipped_steps.name, Status.SKIPPED, 0, None)
        log.error("Run cancelled by the user :(")
        raise

    except Exception as error:
        report(current_step.name, Status.FAILED, 0, None)
        for remaining_step in plan[step_index + 1:]:
            report(remaining_step.name, Status.SKIPPED, 0, None)

        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = files.log_file if files.log_file.is_file() else None
            raise
        raise
