#!/usr/bin/env python3
"""
Testrun with this: python CAT_pack cat -c tests/data/contigs/small_contigs.fa \
                    -d output2/db   -t output2/tax

"""
import shlex
import sys
from decimal import Decimal
from functools import partial
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console, Group
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.text import Text
from typer import Option
from typer_di import Depends, TyperDI

from .pipeline import CatArgs, run_cat, build_plan, Status, Settings, PrepareArgs, run_prepare
from .tools.aligner import DiamondArgs, MMseqsArgs, AlignerName
from .utils.errors import CatError, show_error, InputError
from .utils.logging import init_logging

app = TyperDI()

console = Console(stderr=True)

@app.callback()
def main():
    """Run CAT with progress reporting and collected preflight errors."""



# def show_progress() -> Progress:
#     return Progress(
#         TextColumn("[bold]{task.description:<20}"),
#         BarColumn(bar_width=28, pulse_style="bar.back"),
#         MofNCompleteColumn(separator=" of "),
#         TextColumn("[cyan]{task.fields[status]}"),
#         TimeElapsedColumn(),
#         console=console
# )


class StaticBarColumn(BarColumn):
    def render(self, task):
        bar = super().render(task)
        bar.pulse = False
        return bar


def make_progress(stages):
    progress = Progress(
        TextColumn("[bold]{task.description:<20}"),
        StaticBarColumn(bar_width=28),
        TextColumn("{task.fields[status]}"),
        TimeElapsedColumn(),
        console=console,
    )

    tasks = {
        stage: progress.add_task(
            stage,
            total=1,
            start=False,
            status="[dim]waiting[/dim]",
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
            progress.update(task_id, status=f"[{style}]{status}[/{style}]",
                            refresh=True, total=total, completed=completed)
        case Status.SUPPLIED:
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
        case Status.SKIPPED:
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
        case Status.CANCELLED | Status.FAILED:
            progress.update(task_id, status=f"[{style}]{status}[/{style}]")
            progress.stop_task(task_id)
        case _: # catch all
            progress.update(task_id,status=status,refresh=True)



def system_settings(
        threads: Annotated[
            int,
            Option("--threads", "-n", min=1)
        ] = Settings.threads,
        verbose: Annotated[
            bool,
            Option("--verbose", help="Show aligner stdout."),
        ] = Settings.verbose,
        quiet: Annotated[
            bool,
            Option("--quiet", help="Turns off logging in the terminal."),
        ] = Settings.quiet

):
    return Settings(quiet=quiet, verbose=verbose, threads=threads)

DIAMOND_PATH= Annotated[Path | None, Option(
            "--path-to-diamond", rich_help_panel="DIAMOND",
            help="Path to DIAMOND. Supply if it is not on PATH.")]

@app.command()
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
        db_dir: Annotated[Path,
            Option("--db_dir", help="Directory where CAT/BAT/RAT "
                        "database files will be created", metavar="<DIR>")],
        path_to_diamond: DIAMOND_PATH = None,
        common_prefix: Annotated[str | None, Option(
                "--common-prefix",
                help="Prefix for all files that will be created",
                show_default="<date>_CAT_pack",
            ),] = None,
        cleanup: Annotated["--cleanup"] = None,
        settings: Settings = Depends(system_settings),

):
    args = PrepareArgs(
        db_fasta=db_fasta,
        names=names_dmp,
        nodes=nodes_dmp,
        acc2tax=acc2tax,
        db_dir=db_dir,
        path_to_diamond=path_to_diamond,
        settings=settings,
        common_prefix=common_prefix,
    )
    progress, tasks = make_progress([step.name for step in build_plan(args)])
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




@app.command()
def cat(
        contigs: Annotated[
            Path,
            Option("--contigs", "-c" ,
                   help="Input contig FASTA file", metavar="<file>")
        ],
        database: Annotated[
            Path,
            Option("--database", "-d" ,
                   help="Directory that contains database files",
                   metavar="<directory>")
        ],
        taxonomy: Annotated[
            Path,
            Option("--taxonomy", "-t",
                   help="Directory that contains taxonomy files",
                   metavar="<directory>")
        ],
        range_: Annotated[
            float,
            Option("--range", "-r", min=0.0, max=11,
                   help="r parameter", metavar="<Decimal>"),
        ] = 10.0,
        fraction: Annotated[
            float,
            Option("--fraction", "-f",min=0.0, max=0.99,
                   help="fraction parameter", metavar="<Decimal>"),
        ] = 0.5,
        proteins: Annotated[
            Path | None,
            Option("--proteins_fasta", "-p",
                   help="Predicted proteins fasta file. If supplied, "
                        "the protein prediction step is skipped", metavar="<file>")
        ] = None,
        alignment: Annotated[
            Path | None,
            Option("--alignment_table", "-a",
                   help="Alignment table (in BLAST+6 format). If supplied, "
                    "the alignment step is skipped and classification is "
                    "carried out directly. A predicted proteins fasta file "
                    "should also be supplied with argument --proteins_fasta."
                   , metavar="<file>")
        ] = None,
        output_prefix: Annotated[
            Path,
            Option("--output-prefix", "-o", metavar="<prefix>")
        ] = Path("out.CAT"),
        threads: Annotated[
            int,
            Option("--threads", "-n", min=1)
        ] = Settings.threads,
        top: Annotated[
            int,
            Option("--top", min=0, max=100,
                   help="Hits within range of the best hit written to the alignment file. "
                        "This is not --range."),
        ] = 11,
        tmpdir: Annotated[
            Path | None,
            Option("--tmpdir",
                   help="Location for temporary aligner files."),
        ] = None,
        compress: Annotated[
            bool,
            Option("--compress", help="Compress the alignment output file."),
        ] = False,
        verbose: Annotated[
            bool,
            Option("--verbose", help="Show aligner stdout."),
        ] = False,
        log_file: Annotated[
            Path | None,
            Option("--log-file" , metavar="<file>")
        ] = None,
        debug: Annotated[
            bool,
            Option("--debug", help="Show unexpected-error tracebacks.")
        ] = False,
        aligner: Annotated[
            str,
            Option("--aligner", help="Protein aligner",
                   metavar="<diamond|mmseqs2>", case_sensitive=False)
        ] = "diamond",
        # Seperate Arguments for diamond
        diamond_mode: Annotated[str, Option("--diamond-mode",
            rich_help_panel="DIAMOND",
               help="default, faster, fast, mid-sensitive, sensitive, "
                    "more-sensitive, very-sensitive, ultra-sensitive"
        )] = "default",
        block_size: Annotated[float, Option(
            "--block-size", rich_help_panel="DIAMOND",
            help="DIAMOND block-size. Lower uses less RAM/tmp.",
        )] = 12.0,
        index_chunks: Annotated[int, Option(
            "--index-chunks", min=1, rich_help_panel="DIAMOND",
            help="Set to 4 on low-memory machines.",
        )] = 1,
        no_self_hits: Annotated[bool, Option(
            "--no-self-hits", rich_help_panel="DIAMOND",
            help="Do not report identical self hits by DIAMOND.",
        )] = False,
        path_to_diamond: DIAMOND_PATH = None,
        # Arguments for MMseqs2
        sensitivity: Annotated[float, Option(
            "--sensitivity", min=1.0, max=7.5, rich_help_panel="MMseqs2",
            help="MMseqs2 sensitivity (-s).",
        )] = 5.7,
        split_memory_limit: Annotated[str, Option(
            "--split-memory-limit", rich_help_panel="MMseqs2",
            help="MMseqs2 max memory per split, e.g. 10M, 1G. 0 uses all available memory.",
        )] = "0",
        path_to_mmseqs: Annotated[Path | None, Option(
            "--path-to-mmseqs", rich_help_panel="MMseqs2",
            help="Path to MMseqs2. Supply if it is not on PATH.",
        )] = None,
):

    log = init_logging(debug, quiet=False, log_file=log_file or Path(f"{output_prefix}.log"), console=console)

    log.info("Setting up diamond")
    diamond = DiamondArgs(
        mode=diamond_mode,
        no_self_hits=no_self_hits,
        block_size=block_size,
        index_chunks=index_chunks,
        path_to_diamond=path_to_diamond,
    )

    log.info("Setting up mmseqs")
    mmseqs = MMseqsArgs(
        sensitivity=sensitivity,
        split_memory_limit=split_memory_limit,
        executable=path_to_mmseqs,
    )


    if aligner.lower() == "diamond":
        log.info("Selected aligner is DIAMOND")
        aligner: AlignerName = "diamond"
    elif aligner.lower() in {"mmseqs2", "mmseqs"}:
        log.info("Selected aligner is MMseqs2")
        aligner: AlignerName = "mmseqs2"
    else:
        raise InputError("Aligner must be diamond or mmseqs2")


    log.info("Setting up Arguments for CAT")
    arguments = CatArgs(
        contigs=contigs,
        database=database,
        taxonomy=taxonomy,
        proteins=proteins,
        alignment=alignment,
        range_=Decimal(str(range_)),
        fraction=Decimal(str(fraction)),
        log_file=log_file or Path(f"{output_prefix}.log"),
        output_prefix=output_prefix,
        threads=threads,
        aligner=aligner,
        diamond=diamond,
        mmseqs=mmseqs,
        top=top,
        tmpdir=tmpdir,
        compress=compress,
        verbose=verbose,
    )

    info = Table.grid(padding=(0, 2))
    info.add_column(style="bold cyan")
    info.add_column()

    info.add_row("Contigs", str(arguments.contigs))
    info.add_row("Taxonomy", str(arguments.taxonomy))
    info.add_row("Database", str(arguments.database))
    info.add_row("Aligner", arguments.aligner)
    info.add_row("Parameter r", str(arguments.range_))
    info.add_row("Fraction", str(arguments.fraction))
    info.add_row("Log file", str(arguments.log_file))

    # group the supplied command and parameters together
    content = Group(
        Text("Supplied command", style="bold"),
        Text(f"$ {shlex.join(sys.argv)}", style="cyan"),
        Text(""), info
    )
    log.info(f"Command supplied: $ {shlex.join(sys.argv)}")
    log.info(f"{arguments!r}")

    # print the group in a panel
    console.print(Panel(content, title="[bold]Rarw![/bold]",
                        border_style="blue",), "\n")

    console.print("Preparing for CAT run\n\n")
    log.info("Preparing for CAT run")
    progress, tasks = make_progress([step.name for step in build_plan(arguments)])
    report = partial(update_progress, progress, tasks)
    try:
        with progress:
            outputs = run_cat(arguments, report)
    except KeyboardInterrupt:
        console.print("\nRun cancelled :(")
        raise typer.Exit(code=130)

    except CatError as error:
        show_error(error, console)
        raise typer.Exit(code=1)

    except Exception:
        log.exception("Unexpected error", exc_info=False)
        log.error("Check the run log or use --debug for a full traceback")
        if debug:
            console.print_exception(show_locals=False) # TODO: before release back to False
        raise typer.Exit(code=1)



    else:
        log.info("CAT ran successful!!")
        results = Table(title="CAT completed")
        results.add_column("Result", style="green")
        results.add_column("Location")
        for label, path in outputs.items():
            results.add_row(label, Text(str(path)))
        console.print(results)
