#!/usr/bin/env python3


from typer import Typer, Option
from rich.console import Console




app = Typer()

console = Console()

@app.callback()
def main():
    print("Hi!")

@app.command()
def cat():
    pass