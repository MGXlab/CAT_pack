#!/usr/bin/env python3
import shlex
import sys
from decimal import Decimal
from functools import partial
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console, Group
from rich.panel import Panel
from rich.progress import Progress, ProgressColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.text import Text
from typer import Option

from .cli_options import (
    DiamondPathOption, diamond_options, execution_options, mmseqs_options, with_option_groups,
)
from .config.options import BatOptions, CatOptions, ExecutionOptions, DiamondOptions, MMseqsOptions, PrepareOptions
from .config.validation import get_file_names, make_prefix
from .pipeline import run_annotation, build_plan, Status, run_prepare
from .utils.errors import CatError, show_error
from .utils.logging import init_logging

app = typer.Typer()

console = Console(stderr=True)

@app.callback()
def main():
    """Ah oh"""

class StaticBarColumn(BarColumn):
    def render(self, task):
        bar = super().render(task)
        bar.pulse = False
        return bar


class ProcessingSpeedColumn(ProgressColumn):
    def render(self, task):
        unit = task.fields.get("unit")
        if not unit or not task.started:
            return Text("")
        return Text(f"{task.speed or 0:,.1f} {unit}/s", style="cyan")


def make_progress(stages, unit=None):
    progress = Progress(
        TextColumn("[bold]{task.description:<20}"),
        StaticBarColumn(bar_width=28),
        TextColumn("{task.fields[status]}"),
        TimeElapsedColumn(),
        ProcessingSpeedColumn(),
        console=console,
        refresh_per_second=1
    )

    tasks = {
        stage: progress.add_task(
            stage,
            total=1,
            start=False,
            status="[dim]waiting[/dim]",
            unit=unit if stage == "Classify" else None, # a bit ductaped for now, but for the idea
        )
        for stage in stages
    }

    return progress, tasks


def update_progress(progress, tasks, step, status, completed=0, total=None):
    styles = {
        Status.WAITING: "dim",
        Status.RUNNING: "yellow",
        Status.COMPLETE: "green",
        Status.SUPPLIED: "cyan",
        Status.FAILED: "bold red",
        Status.SKIPPED: "dim",
        Status.CANCELLED: "red",
    }

    style = styles[status]
    task_id = tasks[step]

    # Source https://docs.python.org/3.10/whatsnew/3.10.html#pep-634-structural-pattern-matching
    match status:
        case Status.RUNNING:
            progress.start_task(task_id)
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
            if total is not None:
                progress.update(task_id, total=total, completed=completed)
        case Status.COMPLETE:
            task = next(task for task in progress.tasks if task.id == task_id)
            progress.stop_task(task_id)
            progress.update(task_id, status=f"[{style}]{status}[/{style}]",
                            refresh=True, completed=task.total)
        case Status.SUPPLIED:
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
        case Status.SKIPPED:
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
        case Status.CANCELLED | Status.FAILED:
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
            progress.stop_task(task_id)
        case _: # catch all
            progress.update(task_id,status=status,refresh=True)



@app.command()
@with_option_groups(execution=execution_options)
def prepare(
        db_fasta: Annotated[
            Path, Option("--db_fasta", help="Fasta file containing "
                                            "all sequences.", metavar="<FastaFile>")
        ],
        names_dmp: Annotated[Path, Option("--names",
                                          help="Names.dmp", metavar="<FILE>")],
        nodes_dmp: Annotated[Path, Option("--nodes",
                                          help="Nodes.dmp", metavar="<FILE>")],
        # https://github.com/soedinglab/MMseqs2/blob/master/src/MMseqsBase.cpp#L25-L26
        # MMseqs2 does use this notation of metavar as well, should we look into
        # using this as well?
        acc2tax: Annotated[Path, Option("--acc2tax", help="Accession2taxid.txt file. Can be gzipped.",
                                        metavar="<FILE[.gz]>")],
        # Added --database here but kept db_dir for now too
        db_dir: Annotated[Path,
            Option("--database", "--db_dir", "-d", help="Directory where CAT/BAT/RAT "
                        "database and taxonomy files will be created", metavar="<DIR>")],
        path_to_diamond: DiamondPathOption = PrepareOptions.path_to_diamond,
        common_prefix: Annotated[str | None, Option(
                "--common-prefix",
                help="Prefix for all files that will be created",
                show_default="<date>_CAT_pack",
            ),] = PrepareOptions.common_prefix,
        build_mmseqs2: Annotated[bool, Option("--build-mmseqs2",
                    help="Also create an MMseqs2 sequence database.")] = PrepareOptions.build_mmseqs2,
        path_to_mmseqs: Annotated[Path | None, Option("--path-to-mmseqs",
                    help="Path to MMseqs2, used with --build-mmseqs2.")] = PrepareOptions.path_to_mmseqs,
        *, execution: ExecutionOptions,

):
    args = PrepareOptions(
        db_fasta=db_fasta,
        names=names_dmp,
        nodes=nodes_dmp,
        acc2tax=acc2tax,
        db_dir=db_dir,
        path_to_diamond=path_to_diamond,
        common_prefix=common_prefix,
        build_mmseqs2=build_mmseqs2,
        path_to_mmseqs=path_to_mmseqs,
        threads=execution.threads,
        quiet=execution.quiet,
        verbose=execution.verbose,
        debug=execution.debug,
    )
    args = make_prefix(args)
    progress, tasks = make_progress([step.name for step in build_plan(args)])
    progress.disable = execution.quiet
    report = partial(update_progress, progress, tasks)
    try:
        with progress:
            run_prepare(args, report)
    except KeyboardInterrupt:
        console.print("\nRun cancelled :(")
        raise typer.Exit(code=130)
    except CatError as error:
        show_error(error, console)
        raise typer.Exit(code=1)

    if not execution.quiet:
        files = get_file_names(args)
        table = Table(title="Prepared database")
        table.add_column("File")
        table.add_column("Location")
        for label, path in (
            ("DIAMOND", files.diamond_database),
            ("Protein taxonomy", files.fastaid2LCAtaxid),
            ("Branching taxids", files.taxids_with_multiple_offspring),
            ("Taxonomy names", files.names), ("Taxonomy nodes", files.nodes),
            ("Log", files.log_file),
        ):
            table.add_row(label, str(path))
        if args.build_mmseqs2:
            table.add_row("MMseqs2", str(files.mmseqs2_database))
        console.print(table)
        console.print(f"Use -d {files.data_folder} with CAT_pack7 CAT or BAT.")




@app.command()
@with_option_groups(
    execution=execution_options,
    diamond=diamond_options,
    mmseqs=mmseqs_options,
)
def cat(
        contigs: Annotated[
            Path,
            Option("--contigs", "-c" ,
                   help="Input contig FASTA file", metavar="<file>")
        ],
        database: Annotated[
            Path,
            Option("--database", "-d" ,
                   help="Directory that contains database and taxonomy files",
                   metavar="<directory>")
        ],
        range_: Annotated[
            float,
            Option("--range", "-r", min=0.0, max=100,
                   help="r parameter", metavar="<Decimal>"),
        ] = float(CatOptions.range_),
        fraction: Annotated[
            float,
            Option("--fraction", "-f",min=0.0, max=0.99,
                   help="fraction parameter", metavar="<Decimal>"),
        ] = float(CatOptions.fraction),
        proteins: Annotated[
            Path | None,
            Option("--proteins_fasta", "-p",
                   help="Predicted proteins fasta file. If supplied, "
                        "the protein prediction step is skipped", metavar="<file>")
        ] = CatOptions.proteins,
        alignment: Annotated[
            Path | None,
            Option("--alignment_table", "-a",
                   help="Alignment table (in BLAST+6 format). If supplied, "
                    "the alignment step is skipped and classification is "
                    "carried out directly. A predicted proteins fasta file "
                    "should also be supplied with argument --proteins_fasta."
                   , metavar="<file>")
        ] = CatOptions.alignment,
        output_prefix: Annotated[
            Path,
            Option("--output-prefix", "-o", metavar="<prefix>")
        ] = CatOptions.output_prefix,
        top: Annotated[
            int,
            Option("--top", min=0, max=100,
                   help="Hits within range of the best hit written to the alignment file. "
                        "This is not --range."),
        ] = CatOptions.top,
        tmpdir: Annotated[
            Path | None,
            Option("--tmpdir",
                   help="Location for temporary aligner files."),
        ] = CatOptions.tmpdir,
        compress: Annotated[
            bool,
            Option("--compress", help="Compress the alignment output file."),
        ] = CatOptions.compress,
        log_file: Annotated[
            Path | None,
            Option("--log-file" , metavar="<file>")
        ] = CatOptions.log_file,
        aligner: Annotated[
            str,
            Option("--aligner", help="Protein aligner",
                   metavar="<diamond|mmseqs2>", case_sensitive=False)
        ] = CatOptions.aligner,
        *, execution: ExecutionOptions,
        diamond: DiamondOptions,
        mmseqs: MMseqsOptions,
):

    arguments = CatOptions(
        contigs=contigs,
        database=database,
        proteins=proteins,
        alignment=alignment,
        range_=Decimal(str(range_)),
        fraction=Decimal(str(fraction)),
        log_file=log_file,
        output_prefix=output_prefix,
        aligner=aligner,
        diamond=diamond,
        mmseqs=mmseqs,
        top=top,
        tmpdir=tmpdir,
        compress=compress,
        threads=execution.threads,
        quiet=execution.quiet,
        verbose=execution.verbose, # These two can be place under one name
        debug=execution.debug,     # TODO: merge debug and verbose together
    )

    run_annotation_cli(arguments)

# @bastiaan and @tina, something I found on stackoverflow, for a bit of backwards
# compatibility, we can still let users call CAT_pack bins [OPTIONS] but hide
# it from the --help interface. Yay or Nay?
@app.command("bins", hidden=True)
@app.command()
@with_option_groups(
    execution=execution_options, diamond=diamond_options, mmseqs=mmseqs_options,
)
def bat(
        bins: Annotated[
            Path,
            Option("--bin_fasta", "--bin_folder", "-b",
                   help="Bin fasta file or directory containing bins.", metavar="<FILE|DIR>")
        ],
        database: Annotated[
            Path,
            Option("--database", "-d" ,
                   help="Directory that contains database and taxonomy files",
                   metavar="<directory>")
        ],
        bin_suffix: Annotated[str, Option(
            "--bin_suffix", "--bin-suffix", "-s",
            help="Suffix of bins in bin directory.",
        )] = BatOptions.bin_suffix,
        no_stars: Annotated[bool, Option(
            "--no_stars", "--no-stars",
            help="Suppress marking of suggestive taxonomic assignments.",
        )] = BatOptions.no_stars,
        range_: Annotated[
            float,
            Option("--range", "-r", min=0.0, max=100,
                   help="r parameter", metavar="<Decimal>"),
        ] = float(BatOptions.range_),
        fraction: Annotated[
            float,
            Option("--fraction", "-f",min=0.0, max=0.99,
                   help="f parameter", metavar="<Decimal>"),
        ] = float(BatOptions.fraction),
        proteins: Annotated[
            Path | None,
            Option("--proteins_fasta", "-p",
                   help="Predicted proteins fasta file. If supplied, "
                        "the protein prediction step is skipped.", metavar="<file>")
        ] = BatOptions.proteins,
        alignment: Annotated[
            Path | None,
            Option("--alignment_table", "-a",
                   help="Alignment table (in BLAST+6 format). If supplied, "
                    "the alignment step is skipped and classification is "
                    "carried out directly. A predicted proteins fasta file "
                    "should also be supplied with argument --proteins_fasta."
                   , metavar="<file>")
        ] = BatOptions.alignment,
        output_prefix: Annotated[
            Path,
            Option("--output-prefix", "--out_prefix", "-o", metavar="<prefix>", help="Prefix for output files.")
        ] = BatOptions.output_prefix,
        top: Annotated[
            int,
            Option("--top", min=0, max=100,
                   help="Hits within range of the best hit written to the alignment file. "
                        "This is not --range."),
        ] = BatOptions.top,
        tmpdir: Annotated[
            Path | None,
            Option("--tmpdir",
                   help="Location for temporary aligner files."),
        ] = BatOptions.tmpdir,
        compress: Annotated[
            bool,
            Option("--compress", help="Compress the alignment output file."),
        ] = BatOptions.compress,
        log_file: Annotated[
            Path | None,
            Option("--log-file" , metavar="<file>")
        ] = BatOptions.log_file,
        aligner: Annotated[
            str,
            Option("--aligner", help="Protein aligner",
                   metavar="<diamond|mmseqs2>", case_sensitive=False)
        ] = BatOptions.aligner,
        *, execution: ExecutionOptions,
        diamond: DiamondOptions,
        mmseqs: MMseqsOptions,
):
    """Run Bin Annotation Tool (BAT)."""
    arguments = BatOptions(
        bins=bins,
        bin_suffix=bin_suffix,
        no_stars=no_stars,
        database=database,
        proteins=proteins,
        alignment=alignment,
        range_=Decimal(str(int(range_) if range_ == int(range_) else range_)),
        fraction=Decimal(str(fraction)),
        log_file=log_file,
        output_prefix=output_prefix,
        aligner=aligner,
        diamond=diamond,
        mmseqs=mmseqs,
        top=top,
        tmpdir=tmpdir,
        compress=compress,
        threads=execution.threads,
        quiet=execution.quiet,
        verbose=execution.verbose,
        debug=execution.debug,
    )

    run_annotation_cli(arguments)


def run_annotation_cli(arguments: CatOptions | BatOptions):
    log = init_logging(arguments.debug, quiet=arguments.quiet,
                       console=console)
    log.info("Loaded all arguments")

    info = Table.grid(padding=(0, 2))
    info.add_column(style="bold cyan")
    info.add_column()

    is_bat = isinstance(arguments, BatOptions)
    tool = "BAT" if is_bat else "CAT"
    if is_bat:
        label = "Bin folder" if arguments.bins.is_dir() else "Bin fasta"
        info.add_row(label, str(arguments.bins))
    else:
        info.add_row("Contigs", str(arguments.contigs))
    info.add_row("Database", str(arguments.database))
    info.add_row("Aligner", arguments.aligner)
    info.add_row("Parameter r", str(arguments.range_))
    info.add_row("Fraction", str(arguments.fraction))
    info.add_row("Log file", str(arguments.log_path))

    # group the supplied command and parameters together
    content = Group(
        Text("Supplied command", style="bold"),
        Text(f"$ {shlex.join(sys.argv)}", style="cyan"),
        Text(""), info
    )
    log.info(f"Command supplied: $ {shlex.join(sys.argv)}")
    log.info(f"{arguments!r}")

    # print the group in a panel
    if not arguments.quiet:
        console.print(Panel(content, title="[bold]Rarw![/bold]",
                            border_style="blue",), "\n")


    if not arguments.quiet:
        console.print(f"Preparing for {tool} run\n\n")
    log.info(f"Preparing for {tool} run")
    progress, tasks = make_progress(
        [step.name for step in build_plan(arguments)],
        unit="bins" if is_bat else "contigs",
    )
    progress.disable = arguments.quiet
    report = partial(update_progress, progress, tasks)
    try:
        with progress:
            outputs = run_annotation(arguments, report)
    except KeyboardInterrupt:
        console.print("\nRun cancelled :(")
        raise typer.Exit(code=130)

    except CatError as error:
        show_error(error, console)
        raise typer.Exit(code=1)

    except Exception:
        log.exception("Unexpected error", exc_info=False)
        log.error("Check the run log or use --debug for a full traceback")
        if arguments.debug:
            console.print_exception(show_locals=True) # TODO: before release back to False
        raise typer.Exit(code=1)



    else:
        log.info(f"{tool} ran successfully!")
        if arguments.quiet:
            return

        n_classified, total, entity_type, fraction = outputs["Results"]
        percent = n_classified / total * 100 if total else 0.0
        summary = Text(
            f"{n_classified:,} of {total:,} {entity_type}s "
            f"({percent:.2f}%) have taxonomy assigned.",
            style="bold green",
        )
        results = Table(box=None, padding=(0, 2), expand=True)
        results.add_column("Result", style="green")
        results.add_column("Location", overflow="fold")
        for label, path in outputs.items():
            if isinstance(path, Path):
                results.add_row(Text(label), Text(str(path)))

        content = [summary]
        if fraction < Decimal("0.5"):
            content.append(Text(
                f"Since fraction is set to smaller than 0.5, one {entity_type} "
                f"may have multiple classifications.",
                style="yellow",
            ))
        content.extend([Text("Results are saved at:"), results])
        citation_notes = [step.citation for step in outputs["Steps"] if step.citation]
        if citation_notes:
            content.extend([Text(""), Text("Citation notes", style="bold cyan")])
            content.extend(Text(note) for note in citation_notes)

        console.print(Panel(
            Group(*content),
            title=f"{tool} completed",
            border_style="green",
        ))
