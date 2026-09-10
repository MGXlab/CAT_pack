#!/usr/bin/env python3
from pathlib import Path
from typing import Annotated

import typer

from typer import Typer, Option, Argument
from rich.console import Console

from pipeline import CatArgs



app = Typer()

console = Console()

@app.callback()
def main():
    print("Hi!")


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

):
    print(contigs)
    #arguments = CatArgs(
    #    contigs = contigs)