#!/usr/bin/env python3
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


class CatError(Exception):
    """
    An expected error during a CAT run.

    """
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
        self.log_file: Path | None = None


class InputError(CatError):
    """
    An input error during a CAT run.

    Attributes:
        title: The title of the error.
        hint: The hint of the error.
        path: The path of the error.
        step: The step of the error.

    """
    title = "Check your input"


class ValidationError(CatError):
    title = "Input validation failed"

    def __init__(self, errors: list[CatError]):
        super().__init__(f"Found {len(errors)} problems.")
        self.errors = errors


class ExternalToolError(CatError):
    title = "External tool error"

    def __init__(self, tool: str, message: str, hint: str | None = None):
        self.tool = tool
        super().__init__(f"{tool}: {message}")
        self.hint = hint


class ValErrorCollector:
    """A class that collects all the errors
    Example Usage:
        checks = ValErrorCollector()
        diamond_path = checks.check(check_if_diamond_exists_function, the_path_to_diamond_given by the user)
        Optionally:
            checks.add(SomeError("With other text can be added to the error list"))
        check.finish() # this will raise the VailidationError if any errors where stored

    It can be used nested, that's where the: if isinstance(error, ValidationError)
    comes in to play. That unpacks the "sub" errors of subchecks.

    For example, this can be used with checks for diamond. This way you don't
    need to supply checks, through every function call
    """

    def __init__(self):
        self.errors = []

    def add(self, error):
        if isinstance(error, ValidationError):
            self.errors.extend(error.errors)
        else:
            self.errors.append(error)

    def check(self, function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except CatError as error:
            self.add(error)
            return None
        except OSError as error:
            self.add(InputError(str(error)))
            return None

    def finish(self):
        if self.errors:
            raise ValidationError(self.errors)


def show_error(error: CatError, console: Console) -> None:
    if isinstance(error, ValidationError):
        console.print(
            f"\n[bold red]Input validation failed: "
            f"Found {len(error.errors)} problems.[/bold red]"
        )

        for item in error.errors:
            _show_single_error(item, console)
    else:
        _show_single_error(error, console)


def _show_single_error(error: CatError, console: Console) -> None:
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

    if error.log_file is not None:
        details.add_row("Log", Text(str(error.log_file)))

    console.print(
        Panel(
            details,
            title=Text(error.title, style="bold red"),
            border_style="red",
            expand=True,
            width=90,
            padding=(0, 2),
        )
    )
