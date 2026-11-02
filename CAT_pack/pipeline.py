import logging
from contextlib import ExitStack
from dataclasses import dataclass, replace
from decimal import Decimal
from pathlib import Path
from tarfile import TarError
from tempfile import NamedTemporaryFile
from typing import Any, Callable

from . import tax
from CAT_pack.io.download import (
    DownloadQueueReport, DownloadReport, multi_download,
    check_files, process_files, cleanup_downloads,
)
from CAT_pack.utils.typing import EntityType
from .classification import ClassificationEngine
from .config.options import BatOptions, CatOptions, PrepareOptions, DownloadOptions, RatOptions, AddNamesOptions, \
    SummariseOptions
from .config.settings import BatFiles, BatSettings, CatFiles, CatSettings, DownloadOutputs, AddNamesSettings, \
    SummariseSettings
from .config.validation import check_orfs_match_contigs, get_file_names, get_validated_settings, make_prefix, \
    validate_download, get_download_file_names, make_cat_options, make_bat_options
from .io.parsers import ClassificationParser, optionally_compressed_handle, import_contig_lengths
from .io.writers import ClassificationWriter, write_names, summarise_bins, summarise_contigs
from .prepare import (
    copy_taxonomy, make_diamond_database,
    make_fastaid2LCAtaxid_file, make_mmseqs2_database,
    write_taxids_with_multiple_offspring_file,
)
from .reads import process_bam_file, load_classifications, extract_unclassified_sequences, run_rat_alignment, \
    write_unmapped2classification, calculate_abundances, classify_reads, write_reports
from .tools.aligner import run_aligner
from .tools.bwa_mem import run_bwa_mem
from .tools.pyrodigal import run_protein_prediction
from .utils.check import check_output_prefix, check_outputs
from .utils.errors import CatError
from .utils.locking import lock_outputs
from .utils.logging import Status, Report, file_logging


@dataclass()
class Step:
    """
    A simple dataclass currently acting like a dictionary that stores the
    name of the step and reuse, if reuse is true it will use userprovided data
    """
    name: str
    supplied: bool = False
    citation: str | None = None
    runner: Callable[[Any, Report], Any] | None = None
    result: Any = None
    parent: str | None = None
    display_name: str | None = None
    unit: str | None = None
    total: int | None = 1

    def __str__(self) -> str:
        return self.name


def build_plan(args: CatOptions | BatOptions | RatOptions | PrepareOptions | DownloadOptions | AddNamesOptions, files=None) -> list[Step]:
    if type(args) == PrepareOptions:
        return [
            Step("Input validation"),
            Step("Make DIAMOND database",
                 supplied=files is not None and files.diamond_database.is_file(),
                 runner=make_diamond_database),
            *([Step("Make MMseqs2 database",
                    supplied=files is not None and files.mmseqs2_database.is_file(),
                    runner=make_mmseqs2_database)] if args.build_mmseqs2 else []),
            Step("Make fastaid2LCAtaxid",
                 supplied=files is not None and files.fastaid2LCAtaxid.is_file(),
                 runner=make_fastaid2LCAtaxid_file),
            Step("Make taxids with multiple offspring",
                 supplied=files is not None and files.taxids_with_multiple_offspring.is_file(),
                 runner=write_taxids_with_multiple_offspring_file),
        ]
    if isinstance(args, (CatOptions, BatOptions)):
        return [
            Step("Input validation"),
            Step(
                "Protein prediction",
                supplied=args.proteins is not None,
                citation="Please cite Pyrodigal and Prodigal when using CAT or BAT in your publication.",
                runner=lambda settings, report: run_protein_prediction(
                    settings, report, "pyrodigal"
                ),
            ),
            Step(
                "Alignment",
                supplied=args.alignment is not None,
                runner=lambda settings, report: run_aligner(settings.aligner, report),
                citation="Please cite DIAMOND when using CAT or BAT in your publication."
            ),
            Step("Classify", runner=run_classification,
                 unit="bins" if isinstance(args, BatOptions) else "contigs"),
        ]
    if isinstance(args, RatOptions):
        mode = args.mode.value
        needs_alignment = ("r" in mode and args.unmapped2classification is None)
        needs_sequences = (needs_alignment and args.alignment_unmapped is None)

        return [
            Step("Input validation"),
            Step("Read mapping", supplied=bool(args.bam_file), runner=run_bwa_mem),
            *(annotation_steps(make_cat_options(args), "cat", "Contigs",
                               args.contigs2classification is not None) if "c" in mode else []),
            *(annotation_steps(make_bat_options(args), "bat", "Bins",
                               args.bin2classification is not None) if "m" in mode else []),
            (proc_maps := Step("Process mappings", runner=process_bam_file)),
            (classif := Step("Load classifications", total=4,
                             runner=lambda settings, report:
            load_classifications(settings, report, proc_maps.result))),
            *([Step("Extract unclassified sequences",
                    runner=lambda settings, report: extract_unclassified_sequences(
            settings, report, classif.result))] if needs_sequences else []),
            *([Step("Align unclassified sequences",supplied=args.alignment_unmapped is not None, runner=run_rat_alignment)] if needs_alignment else []),
            *([Step("Classify unclassified sequences",supplied=args.unmapped2classification is not None,  total=4,
                runner=lambda settings, report: write_unmapped2classification(settings, report, classif.result)
                    )] if "r" in mode else []),
            Step(
                "Calculate abundances",
                total=4,
                runner=lambda settings, report:
                calculate_abundances(
                    settings, report, classif.result,
                ),
            ),
            *([Step(
                "Assign read classifications",
                runner=lambda settings, report:
                classify_reads(
                    settings, report, classif.result,
                ),
            )] if args.read2classification else []),
            Step(
                "Write reports",
                total=3,
                runner=lambda settings, report:
                write_reports(
                    settings, report, classif.result,
                ),
            ),
        ]
    if isinstance(args, DownloadOptions):
        return [
            Step("Input validation"),
            Step("Downloading"),
            Step("Check files"),
            Step("Process files"),
            Step("Cleanup", supplied=not args.cleanup),
        ]
    if isinstance(args, AddNamesOptions):
        return [
            Step("Input validation"),
            Step("Add names", total=3, runner=append_names),
        ]
    if isinstance(args, SummariseOptions):
        return [
            Step("Input validation"),
            Step("Summarise", total=2 if args.contigs_fasta else 1, runner=summarise)
        ]
    raise CatError("I haven't figured out how to build that specific plan")


def run_classification(settings: CatSettings | BatSettings, report: Report) -> tuple[int, int, EntityType, Decimal]:
    """Load inputs, classify all entities, then write their reports."""
    files: CatFiles | BatFiles = settings.files
    log = logging.getLogger("CAT_pack")
    inputs = ClassificationParser(
        files, settings.range_, parser_settings=settings.alignment_parsing, workers=settings.threads
    ).parse()
    check_orfs_match_contigs(
        inputs.contig_names, inputs.contig2ORFs, files.proteins_fasta
    )
    engine = ClassificationEngine(
        taxonomy=inputs.taxonomy,
        fastaid2taxid=inputs.fastaid2taxid,
        fraction=settings.fraction,
        threads=settings.threads,
    )

    total = len(inputs.entity2ORFs)
    report("Classify", Status.RUNNING, 0, total)
    result = engine.classify(
        inputs.entity_type, inputs.entity2ORFs, inputs.alignment
    )

    log.info(
        "Writing classifications to %s and %s.",
        files.report, files.orf_report,
    )
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
            no_stars=settings.no_stars,
        )
        for entity in range(len(result)):
            writer.write(result, entity)
            report("Classify", Status.RUNNING, entity + 1, total)

    n_classified = result.n_assigned
    classification_summary(
        n_classified, total, inputs.entity_type, settings.fraction
    )
    return n_classified, total, inputs.entity_type, settings.fraction


def classification_summary(n_classified: int, total: int, entity_type: str, fraction: Decimal) -> None:
    log = logging.getLogger("CAT_pack")
    tool = "CAT"
    if entity_type == EntityType.BIN:
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

def run_download(
    args: DownloadOptions, report: Report,
    download_report: DownloadReport,
    queue_report: DownloadQueueReport
) -> DownloadOutputs:
    files = get_download_file_names(args)
    plan = build_plan(args)
    step_index = 0
    current_step = plan[step_index]
    log = logging.getLogger("CAT_pack")
    locks = ExitStack()
    settings = None
    try:
        files.data_folder.mkdir(parents=True, exist_ok=True)
        locks.enter_context(lock_outputs([files.data_folder / ".CAT_pack.download.lock"]))
        for step_index, current_step in enumerate(plan):
            # A bit of a misuse of the term supplied. this is for if cleanup isn't enabled
            if current_step.supplied:
                report(current_step.name, Status.SKIPPED, 0, None)
                continue
            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")
            if step_index == 0:
                settings = validate_download(args, files)
                if not settings.no_log:
                    locks.enter_context(file_logging(log, files.log_file, settings.debug))
            elif current_step.name == "Downloading":
                multi_download(settings, report, download_report, queue_report)
            elif current_step.name == "Check files":
                check_files(files, report)
            elif current_step.name == "Process files":
                process_files(files, report)
            elif current_step.name == "Cleanup":
                cleanup_downloads(files)

            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")

        log.info(f"Download complete. {settings.db_type.upper()} inputs are ready for CAT_pack7 prepare.")
        return files

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
        log.error("Download failed during %s: %s", current_step.name, error)

        if isinstance(error, (OSError, UnicodeError, EOFError, TarError)):
            error = CatError(str(error), path=files.data_folder)
        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = files.log_file if not args.no_log and files.log_file.is_file() else None
            raise error
        raise
    finally:
        locks.close()

def annotation_steps(options, settings_name: str, label: str, supplied: bool) -> list[Step]:
    if supplied:
        return [Step(f"Classify {label.lower()}", supplied=True)]

    def adapt(step):
        name = f"{label}: {step.name}"

        def runner(settings, report):
            def forward_progress(_, status, completed, total):
                report(name, status, completed, total)

            return step.runner(getattr(settings, settings_name), forward_progress)

        return replace(step, name=name, runner=runner,
                       parent=f"Classify {label.lower()}", display_name=step.name)

    return [adapt(step) for step in build_plan(options)[1:]]


def run(args: CatOptions | BatOptions | RatOptions | PrepareOptions | AddNamesOptions | SummariseOptions, report: Report):
    log = logging.getLogger("CAT_pack")

    is_prepare = isinstance(args, PrepareOptions)
    is_add_names = isinstance(args, AddNamesOptions)
    is_summarise = isinstance(args, SummariseOptions)
    if is_prepare: args = make_prefix(args)
    log_path = get_file_names(args).log_file if is_prepare else args.log_path

    plan = build_plan(args)
    completed_steps: list[Step] = []
    step_index = 0
    current_step = plan[step_index]
    locks = ExitStack()

    try:

        log.info("Doing some pre-flight checks first.")
        report(current_step.name, Status.RUNNING, 0, 1)
        log.info(f"Starting: {current_step.name}")
        if is_prepare:
            args.db_dir.mkdir(parents=True, exist_ok=True)
            locks.enter_context(lock_outputs([args.db_dir / ".prepare.lock"]))
        else:
            if is_add_names or is_summarise:
                prefixes = [args.output_file]
            else:
                prefixes = [args.output_prefix]
            if isinstance(args, RatOptions):
                if "c" in args.mode.value and args.contigs2classification is None:
                    prefixes.append(make_cat_options(args).output_prefix)
                if "m" in args.mode.value and args.bin2classification is None:
                    prefixes.append(make_bat_options(args).output_prefix)
            outputs = [*prefixes, log_path]
            for path in outputs:
                check_output_prefix(path)
            locks.enter_context(lock_outputs(path.with_name(f".{path.name}.lock") for path in outputs))

        settings = get_validated_settings(args)

        if is_prepare:
            locks.enter_context(lock_outputs([log_path.with_name(f".{log_path.name}.lock")]))
        locks.enter_context(file_logging(log, log_path, args.debug))
        if is_prepare: copy_taxonomy(settings)

        plan = build_plan(args, settings.files)
        report(current_step.name, Status.COMPLETE, 1, 1)
        log.info(f"Completed: {current_step.name}")

        for step_index, current_step in enumerate(plan[1:], start=1):
            if current_step.supplied:
                report(current_step.name, Status.SUPPLIED, 1, 1)
                log.info(f"Supplied file: {current_step.name}")
                continue
            report(current_step.name, Status.RUNNING, 0, None)
            log.info(f"Starting: {current_step.name}")
            if current_step.runner is None:
                raise CatError(f"{current_step.name} is not implemented yet.")
            current_step.result = current_step.runner(settings, report)
            report(current_step.name, Status.COMPLETE, 1, 1)
            log.info(f"Completed: {current_step.name}")
            completed_steps.append(current_step)
        if is_prepare:
            log.info("Preparation complete. Use -d %s with CAT_pack7 CAT or BAT.", settings.files.data_folder)
            return 0
        if isinstance(args, (SummariseOptions, AddNamesOptions)):
            return {
                "Results": completed_steps[-1].result,
                "Steps": completed_steps,
                "Named output": settings.files.output_file,
                "Log": settings.log_file,
            }
        if isinstance(args, RatOptions):
            return {
                "Steps": completed_steps,
                "Contig abundance": settings.files.contig_abundance,
                "Complete abundance": settings.files.complete_abundance,
                **(
                    {"Bin abundance": settings.files.bin_abundance}
                    if settings.files.bin_abundance is not None else {}
                ),
                **(
                    {"Read classifications": settings.files.read_classification}
                    if settings.files.read_classification is not None else {}
                ),
                "Log": settings.log_file,
                "Results": completed_steps[-1].result,
            }
        return {
            "Results": completed_steps[-1].result,
            "Steps": completed_steps,
            "Classifications": settings.files.report,
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

        if isinstance(error, (OSError, UnicodeError, EOFError)):
            error = CatError(str(error), path=args.db_dir if is_prepare else args.output_file if is_add_names or is_summarise else args.output_prefix)
        if isinstance(error, CatError):
            error.step = current_step.name
            error.log_file = log_path
            raise error
        # TODO: REMOVE incomplete files
        # Catch all other exceptions
        raise
    finally:
        locks.close()


def append_names(settings: AddNamesSettings, report: Report) -> Path:
    log = logging.getLogger("CAT_pack")
    files = settings.files

    taxid2rank = tax.read_nodes(files.nodes, ranks=True)
    report.progress_next()

    taxid2name = tax.import_names(files.names)
    report.progress_next()

    log.info("Appending names...")
    temporary = None

    try:
        with NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=files.output_file.parent,
            prefix=f".{files.output_file.name}.",
            suffix=".tmp",
            delete=False,
        ) as outf1:
            temporary = Path(outf1.name)

            with optionally_compressed_handle(files.input_file) as f1:
                write_names(
                    f1,
                    outf1.file,
                    taxid2rank,
                    taxid2name,
                    input_file=files.input_file,
                    only_official=settings.only_official,
                    exclude_scores=settings.exclude_scores,
                )

        if not settings.force:
            check_outputs([files.output_file])
        temporary.replace(files.output_file)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

    report.progress_next()
    log.info(f"Names written to {settings.files.output_file.name}")
    return files.output_file

def summarise(settings: SummariseSettings, report: Report) -> Path:
    log = logging.getLogger("CAT_pack")
    files = settings.files
    contig2length = {}

    if files.contigs_fasta is not None:
        contig2length = import_contig_lengths(files.contigs_fasta)
        report.progress_next()

    log.info("Summarising...")
    temporary = None

    try:
        with NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=files.output_file.parent,
            prefix=f".{files.output_file.name}.",
            suffix=".tmp",
            delete=False,
        ) as outf1:
            temporary = Path(outf1.name)

            with optionally_compressed_handle(files.input_file) as f1:
                if files.contigs_fasta is None:
                    summarise_bins(
                        f1, outf1.file,
                        input_file=files.input_file,
                    )
                else:
                    summarise_contigs(
                        f1,
                        outf1.file,
                        contig2length,
                        input_file=files.input_file,
                        contigs_fasta=files.contigs_fasta,
                    )

        if not settings.force:
            check_outputs([files.output_file])
        temporary.replace(files.output_file)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

    report.progress_next()
    return files.output_file
