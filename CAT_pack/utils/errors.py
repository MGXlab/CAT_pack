#!/usr/bin/env python3
from pathlib import Path
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

class CatError(Exception):
    """An expected error during a CAT run."""
    title = "CAT could not finish"
    exit_code = 1

    def __init__(self, message: str, *,
            hint: str | None = None,
            path: Path | None = None
    ):
        super().__init__(message)

        self.hint = hint
        self.path = path
        self.step: str | None = None
        #self.log_file: Path | None = None


class InputError(CatError):
    title = "Check your input"


class ExternalToolError(CatError):
    def __init__(self, tool: str, message: str):
        self.tool = tool
        super().__init__(f"{tool}: {message}")



def show_error(error: CatError, console: Console) -> None:
    details = Table.grid(padding=(0, 3))
    details.add_column(style="bold", no_wrap=True)
    details.add_column()

    details.add_row("Error message:", Text(str(error)))

    if error.step is not None:
        details.add_row("Step", Text(error.step))

    if error.path is not None:
        details.add_row("Path", Text(str(error.path)))

    if error.hint is not None:
        details.add_row("Try this", Text(error.hint))

    console.print(
        Panel(
            details,
            title=Text(error.title, style="bold red"),
            border_style="red",
            expand=False,
            padding=(0, 2),
        )
    )