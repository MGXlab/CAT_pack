#!/usr/bin/env python3
"""
Testrun with this: python CAT_pack cat -c tests/data/contigs/small_contigs.fa   -d output2/db   -t output2/tax   -n 4   -o CAT_run_new   --verbose --force

"""



from decimal import Decimal
from functools import partial
from pathlib import Path
from typing import Annotated

import typer

from typer import Typer, Option, Argument
from rich.console import Console
from rich.progress import Progress, BarColumn, MofNCompleteColumn, TextColumn, TimeElapsedColumn

from pipeline import CatArgs, run_cat, stages




app = Typer()

console = Console(stderr=True)

@app.callback()
def main():
    print("Hi!")



# def show_progress() -> Progress:
#     return Progress(
#         TextColumn("[bold]{task.description:<20}"),
#         BarColumn(bar_width=28, pulse_style="bar.back"),
#         MofNCompleteColumn(separator=" of "),
#         TextColumn("[cyan]{task.fields[status]}"),
#         TimeElapsedColumn(),
#         console=console
# )


def make_progress(stages):
    progress = Progress(
        TextColumn("[bold]{task.description:<20}"),
        BarColumn(bar_width=28, pulse_style="bar.back"),
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


def update_progress(progress, tasks, stage, status, completed=0, total=None):
    task_id = tasks[stage]

    if status == "running":
        progress.start_task(task_id)

        if total is not None:
            progress.update(task_id, total=total, completed=completed)

    elif status == "complete":
        progress.update(task_id, total=total, completed=completed)

    elif status in {"failed", "cancelled"}:
        progress.stop_task(task_id)

    progress.update(task_id,status={status},refresh=True)




# @bastiaan there is a destinction between typer.Option and typer.Arguments
# Arguments are stricly bound to input order of arguments, and do not allow for aliases
# But they are by default required
# Options on the other hand must be set by a argument name and allow for aliases
# However they are by default NOT required. But can be set to be required
# What do you think is best? For now I will go for Options and we can always re-evaluate
@app.command()
def cat(
        contigs: Annotated[
            Path,
            Option("--contigs", "-c" ,help="Input contig FASTA file" )
        ],
        database: Annotated[
            Path,
            Option("--database", "-d" ,help="Directory that contains database files" )
        ],
        taxonomy: Annotated[
            Path,
            Option("--taxonomy", "-t", help="Directory that contains taxonomy files" )
        ],
        range_: Annotated[
            float,
            Option("--range", "-r", min=0.0, max=100, help="r parameter"),
        ] = 10.0,
        proteins: Annotated[
            Path,
            Option("--proteins_fatsa", "-p",
                   help="Predicted proteins fasta file. If supplied, "
                        "the protein prediction step is skipped"),
        ] = None
):
    # notes: Decimal is not supported by typer (look into that)
    # Print is only for my own debugging for now
    #print(contigs, database, taxonomy)

    arguments = CatArgs(
        contigs=contigs,
        database=database,
        taxonomy=taxonomy,
        proteins=proteins,
        range_=Decimal(str(range_)),
        output_prefix=Path("./out.CAT"),
    )

    #console.print(arguments)
    console.print("Ready for takeoff")
    progress, tasks = make_progress(stages)
    report = partial(update_progress, progress, tasks)
    try:
        with progress:
            outputs = run_cat(arguments, report)
    except KeyboardInterrupt:
        console.print("\nRun cancelled :(")
        raise typer.Exit(code=130)

    except Exception as error:
        console.print(error)
        raise typer.Exit(code=1)



    else:
        console.print(outputs)