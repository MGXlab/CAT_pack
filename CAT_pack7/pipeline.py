import logging
from contextlib import ExitStack
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from . import tax
from .classification import ClassificationEngine
from .config.options import BatOptions, CatOptions, PrepareOptions
from .config.settings import BatFiles, BatSettings, CatFiles, CatSettings, DiamondSettings
from .config.validation import check_orfs_match_contigs, get_file_names, get_validated_settings, make_prefix, \
    validate_prepare
from .io.parsers import ClassificationParser
from .io.writers import ClassificationWriter
from .prepare import (
    copy_taxonomy, find_offspring, make_diamond_database,
    make_fastaid2LCAtaxid_file, make_mmseqs2_database,
    write_taxids_with_multiple_offspring_file,
)
from .tools.aligner import run_aligner
from .tools.pyrodigal import run_protein_prediction
from .utils.errors import CatError
from .utils.locking import lock_outputs
from .utils.logging import Status, Report, file_logging


@dataclass()
class Step:
    """
    A simple dataclass currently acting like a dictionairy that stores the
    name of the step and reuse, if reuse is true it will use userprovided data
    """

    name: str
    supplied: bool = False
    citation: str | None = None


    def __str__(self) -> str:
        return self.name


def build_plan(args: CatOptions | BatOptions | PrepareOptions, files=None) -> list[Step]:
    if type(args) == PrepareOptions:
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
    if isinstance(args, (CatOptions, BatOptions)):
        return [
            Step("Input validation"),
            Step(
                "Protein prediction", supplied=args.proteins is not None,
                citation="Please cite Pyrodigal and Prodigal when using CAT or BAT in your publication.",
            ),
            Step("Alignment", supplied=args.alignment is not None),
            Step("Classify"),
        ]
    raise CatError("I haven't figured out how to build that specific plan")

def run_annotation(
    args: CatOptions | BatOptions, report: Report,
) -> dict[str, Path | tuple[int, int, str, Decimal] | list[Step]]:
    """Bin/Contig annotation tool (BAT/CAT) run"""
    is_bat = isinstance(args, BatOptions)

    plan = build_plan(args)
    completed_steps: list[Step] = []
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")
    locks = ExitStack()

    try:
        locks.enter_context(lock_outputs([args.output_prefix, args.log_path]))
        if args.alignment is None:
            tmpdir = args.tmpdir if args.tmpdir is not None else args.output_prefix.parent / "tmp"
            tmpdir.parent.mkdir(parents=True, exist_ok=True)
            locks.enter_context(lock_outputs([tmpdir]))

        locks.enter_context(file_logging(log, args.log_path, args.debug))

        planned_steps = ", ".join(step.name for step in plan[1:] if not step.supplied)
        log.info(f"{'BAT' if is_bat else 'CAT'} is running. Planned steps: {planned_steps}.")
        log.info("Doing some pre-flight checks first.")
        report(current_step.name, Status.RUNNING, 0, 1)
        log.info(f"Starting: {current_step.name}")
        settings = get_validated_settings(args)
        report(current_step.name, Status.COMPLETE, 1, 1)
        log.info(f"Completed: {current_step.name}")
        completed_steps.append(current_step)
        if is_bat:
            log.info("Ready to fly!\n\n-----------------\n")

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
                if isinstance(settings.aligner, DiamondSettings):
                    current_step.citation = (
                        "Please cite DIAMOND when using CAT or BAT in your publication."
                    )
            elif current_step.name == "Classify":
                results = run_classification(settings, settings.files, report)
            else:
                log.error(f"Can't spell that well, my misspelling:"
                          f"{current_step.name}")

            # always complete the step with 1 of 1 so 100% completion,
            # errors will have caused the code to stop before this if there is
            # an exception.
            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")
            completed_steps.append(current_step)

        return {
            "Results": results,
            "Steps": completed_steps,
            f"{'Bin' if is_bat else 'Contig'} Classifications": settings.files.report,
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

        if isinstance(error, OSError):
            error = CatError(str(error), path=args.output_prefix)
        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = args.log_path if args.log_path.is_file() else None
            raise error

        # Catch all other exceptions
        raise
    finally:
        locks.close()


def run_classification(
    settings: CatSettings | BatSettings,
    files: CatFiles | BatFiles,
    report: Report,
) -> tuple[int, int, str, Decimal]:
    """Load and validate inputs, then classify and write the results."""
    log = logging.getLogger("CAT_pack")
    inputs = ClassificationParser(files, settings.range_).parse()
    check_orfs_match_contigs(inputs.contig_names, inputs.contig2ORFs, files.proteins_fasta)
    engine = ClassificationEngine(
        taxid2parent=inputs.taxid2parent,
        fastaid2taxid=inputs.fastaid2taxid,
        fraction=settings.fraction,
    )

    if isinstance(files, BatFiles):
        tool = "BAT"
        action = "flying"
    else:
        tool = "CAT"
        action = "spinning"
    log.info(f"{tool} is {action}! Files {files.report} and {files.orf_report} are created.")
    n_classified = 0
    total = len(inputs.entity2ORFs)
    report("Classify", Status.RUNNING, 0, total)
    with (
        files.report.open("w", encoding="utf-8") as classification_out,
        files.orf_report.open("w", encoding="utf-8") as orf_out,
    ):
        writer = ClassificationWriter(
            classification_out,
            orf_out,
            entity_type=inputs.entity_type,
            fraction=settings.fraction,
            range_=settings.range_,
            branches=inputs.branches,
            no_stars=getattr(settings, "no_stars", False),
        )
        for index, entity_id in enumerate(sorted(inputs.entity2ORFs), start=1):
            result = engine.classify_group(
                entity_id=entity_id,
                orf_ids=inputs.entity2ORFs[entity_id],
                orf2hits=inputs.alignment.orf2hits,
            )
            writer.write(result)
            if result.assigned:
                n_classified += 1
            report("Classify", Status.RUNNING, index, total)
    classification_summary(n_classified, total, inputs.entity_type, settings.fraction)
    return n_classified, total, inputs.entity_type, settings.fraction


def classification_summary(
    n_classified: int,
    total: int,
    entity_type: str,
    fraction: Decimal,
) -> None:
    log = logging.getLogger("CAT_pack")
    tool = "CAT"
    if entity_type == "bin":
        tool = "BAT"
    percent = 0
    if total > 0:
        percent = n_classified / total * 100
    log.info(
        f"{tool} is done! {n_classified:,d}/{total:,d} {entity_type}s "
        f"({percent:.2f}%) have taxonomy assigned."
    )
    if fraction < Decimal("0.5"):
        log.warning(
            f"since f is set to smaller than 0.5, one {entity_type} may have multiple classifications."
        )


def run_prepare(args: PrepareOptions, report: Report):
    args = make_prefix(args)
    files = get_file_names(args)
    plan = build_plan(args, files)
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")
    locks = ExitStack()
    taxid2parent = None

    try:
        files.data_folder.parent.mkdir(parents=True, exist_ok=True)
        locks.enter_context(lock_outputs([files.data_folder]))
        plan = build_plan(args, files)
        for step_index, current_step in enumerate(plan):
            if current_step.name == "Make MMseqs2 database" and not args.build_mmseqs2:
                report(current_step.name, Status.SKIPPED, 0, None)
                log.info("MMseqs2 database creation was not requested.")
                continue
            if current_step.supplied:
                report(current_step.name, Status.SUPPLIED, 1, 1)
                log.info(f"Already exists, skipped making of: {current_step.name}")
                continue

            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")
            if step_index == 0:
                settings = validate_prepare(args, files)
                files.data_folder.mkdir(parents=True, exist_ok=True)
                locks.enter_context(file_logging(log, files.log_file, args.debug))
                copy_taxonomy(settings)
            elif current_step.name == "Make DIAMOND database":
                make_diamond_database(settings)
            elif current_step.name == "Make MMseqs2 database":
                make_mmseqs2_database(settings)
            else:
                if taxid2parent is None:
                    taxid2parent, _ = tax.import_nodes(files.nodes)
                if current_step.name == "Make fastaid2LCAtaxid":
                    make_fastaid2LCAtaxid_file(
                        files.fastaid2LCAtaxid, settings.db_fasta, settings.acc2tax,
                        taxid2parent, report,
                    )
                elif current_step.name == "Make taxids with multiple offspring":
                    taxid2offspring = find_offspring(files.fastaid2LCAtaxid, taxid2parent, report)
                    write_taxids_with_multiple_offspring_file(
                        files.taxids_with_multiple_offspring, taxid2offspring,
                    )
            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")

        log.info("Preparation complete. Use -d %s with CAT_pack7 CAT or BAT.", files.data_folder)
        return 0

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
        log.error("Preparation failed during %s: %s", current_step.name, error)

        if isinstance(error, (OSError, UnicodeError, EOFError)):
            error = CatError(str(error), path=files.data_folder)
        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = files.log_file if files.log_file.is_file() else None
            raise error
        raise
    finally:
        locks.close()
