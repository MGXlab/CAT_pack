import logging
from dataclasses import dataclass
from datetime import datetime
from decimal import Decimal
from pathlib import Path

from .classification import contig_classification
from .tools.aligner import run_aligner, DiamondArgs, MMseqsArgs, AlignerName
from .tools.pyrodigal import run_protein_prediction
from .utils.errors import CatError
from .utils.logging import Status, Report
from .validation import validate_args, validate_aligner_args


# Defaults
# Parent

# CatDefault
# Child


@dataclass(kw_only=True)
class Defaults:
    quiet: bool = False
    verbose: bool = False
    threads: int = 1
    db: Path = None
    path_to_diamond: Path = None
    path_to_mmseqs2: Path = None
    common_prefix: Path | Path = f"{datetime.now():%Y-%m-%d}_CAT_pack"
    output_prefix: Path | None = None
    log_file: Path | Path = f"{datetime.now():%Y-%m-%d}_CAT_pack.log"



@dataclass(kw_only=True)
class PrepareDefaults(Defaults):
    db_fasta: Path
    names: Path
    nodes: Path
    acc2tax: Path
    db_dir: Path
    path_to_diamond: Path | None
    threads: int = 1


@dataclass(kw_only=True)
class CatDefaults(Defaults):
    contigs: Path
    database: Path
    taxonomy: Path
    proteins: Path | None = None
    alignment: Path | None = None
    range_: Decimal = 10.0
    fraction: Decimal = 0.5
    diamond: DiamondArgs
    mmseqs: MMseqsArgs
    aligner: AlignerName = "diamond"
    top: int = 11
    tmpdir: Path | None = None
    compress: bool = False



@dataclass()
class Settings: #@Bastiaan Know a better name?
    quiet: bool = False
    verbose: bool = False
    threads: int = 1


@dataclass(frozen=True)
class CatArgs:
    contigs: Path
    database: Path
    taxonomy: Path
    proteins: Path | None
    alignment: Path | None
    range_: Decimal
    fraction: Decimal
    log_file: Path
    output_prefix: Path
    diamond: DiamondArgs
    mmseqs: MMseqsArgs
    aligner: AlignerName = "diamond"
    threads: int = 1
    top: int = 11
    tmpdir: Path | None = None
    compress: bool = False
    verbose: bool = False

@dataclass(frozen=True)
class PrepareArgs:
    def __init_subclass__(cls, **kwargs):
        pass


    db_fasta: Path
    names: Path
    nodes: Path
    acc2tax: Path
    db_dir: Path
    path_to_diamond: Path | None
    defaults: Defaults


@dataclass(frozen=True)
class PrepareOutputs:
    prefix: str
    db_folder: Path
    tax_folder: Path
    log_file: Path
    diamond_database: Path
    mmseqs2_database: Path
    fastaid2LCAtaxid: Path
    taxids_with_multiple_offspring: Path


def expand_prepare(args: PrepareArgs) -> PrepareOutputs:
    """Mini expand_arguments van shared.py"""
    prefix = args.defaults.common_prefix or f"{datetime.now():%Y-%m-%d}_CAT_pack"
    db_folder = args.db_dir / "db"
    tax_folder = args.db_dir / "tax"
    return PrepareOutputs(
        prefix=prefix,
        db_folder=db_folder, # on folder
        tax_folder=tax_folder,
        log_file=args.db_dir / f"{prefix}.log",
        diamond_database=db_folder / f"{prefix}.dmnd",
        mmseqs2_database=db_folder / f"{prefix}.mmseqs2",
        fastaid2LCAtaxid=db_folder / f"{prefix}.fastaid2LCAtaxid",
        taxids_with_multiple_offspring=db_folder / f"{prefix}.taxids_with_multiple_offspring",
    )

@dataclass(frozen=True)
class Step:
    """
    A simple dataclass currently acting like a dictionairy that stores the
    name of the step and reuse, if reuse is true it will use userprovided data
    """

    name: str
    supplied: bool = False

    def __str__(self) -> str:
        return self.name


def build_plan(args: CatArgs | PrepareArgs) -> list[Step]:
    if type(args) == PrepareArgs:
        files = expand_prepare(args)
        return [
            Step("Input validation"),
            Step("Make DIAMOND database", supplied=files.diamond_database.is_file()),
            Step("Make MMseqs2 database", supplied=files.mmseqs2_database.is_file()),
            Step("Make fastaid2LCAtaxid", supplied=files.fastaid2LCAtaxid.is_file()),
            Step("Make taxids with multiple offspring",
                 supplied=files.taxids_with_multiple_offspring.is_file()),
        ]
    if type(args) == CatArgs:
        return [
            Step("Input validation"),
            Step("Protein prediction", supplied=args.proteins is not None),
            Step("Alignment", supplied=args.alignment is not None),
            Step("Classify"),
        ]
    raise CatError("I haven't figured out how to build that specific plan")


def run_cat(args: CatArgs, report: Report) -> dict[str, Path]:
    """Contig annotation tool (CAT) run"""

    plan = build_plan(args)
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")


    try:
        report(current_step.name, Status.RUNNING, 0, 1)
        log.info(f"Starting: {current_step.name}")
        files = validate_args(args) # TODO: still need to revamp this
        aligner_args = validate_aligner_args(args, files)
        report(current_step.name, Status.COMPLETE, 1, 1)
        log.info(f"Completed: {current_step.name}")

        for step_index, current_step in enumerate(plan[1:], start=1):
            if current_step.supplied:
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





def run_prepare(args: PrepareArgs, report: Report):
    plan = build_plan(args)
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")
    files = expand_prepare(args)

    try:
        for step_index, current_step in enumerate(plan):
            if current_step.supplied:
                report(current_step.name, Status.SUPPLIED, 1, 1)
                log.info(f"Already exists, skipped making of: {current_step.name}")
                continue

            if step_index == 0:
                # TODO: validate prepare inputs (fasta, names, nodes, acc2tax)
                pass

            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")
            # TODO: Port over the actual steps
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
            error.log_file = files.log_file
            raise
        raise