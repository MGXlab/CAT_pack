#!/usr/bin/env python3
"""
Testrun with this: python CAT_pack cat -c tests/data/contigs/small_contigs.fa   -d output2/db   -t output2/tax   -n 4   -o CAT_run_new   --verbose --force

"""



from decimal import Decimal
from pathlib import Path
from typing import Annotated

import typer

from typer import Typer, Option, Argument
from rich.console import Console
from rich.progress import Progress, BarColumn, MofNCompleteColumn, TextColumn, TimeElapsedColumn

from pipeline import CatArgs




app = Typer()

console = Console(stderr=True)

@app.callback()
def main():
    print("Hi!")


# @bastiaan there is a destinction between typer.Option and typer.Arguments
# Arguments are stricly bound to input order of arguments, and do not allow for aliases
# But they are by default required
# Options on the other hand must be set by a argument name and allow for aliases
# However they are by default NOT required. But can be set to be required
# What do you think is best? For now I will go for Options and we can always re-evaluate
def show_progress() -> Progress:
    return Progress(
        TextColumn("[bold]{task.description:<20}"),
        BarColumn(bar_width=28, pulse_style="bar.back"),
        MofNCompleteColumn(separator=" of "),
        TextColumn("[cyan]{task.fields[status]}"),
        TimeElapsedColumn(),
        console=console
)


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
):
    # notes: Decimal is not supported by typer (look into that)
    # Print is only for my own debugging for now
    print(contigs, database, taxonomy)

    arguments = CatArgs(
        contigs=contigs,
        database=database,
        taxonomy=taxonomy,
        _range=Decimal(str(range_)),
        output_prefix=Path("./out.CAT"),
    )

    console.print(arguments)
    console.print("Ready for takeoff")
    progress = show_progress()
    try:
        outputs = run_cat(
            arguments,
            on_event=progress.handle_event,
        )
    except KeyboardInterrupt:
        console.print("\nRun cancelled :(")
        raise typer.Exit(code=130)

    except Exception as error:
        console.print(error)
        raise typer.Exit(code=1)



    else:
        console.print(outputs)