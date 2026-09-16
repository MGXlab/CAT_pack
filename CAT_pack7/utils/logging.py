import logging
from enum import Enum
from pathlib import Path
from typing import Protocol

from rich.console import Console
from rich.logging import RichHandler


class Status(Enum):
    WAITING = "Waiting"
    RUNNING = "Running"
    COMPLETE = "Completed"
    SUPPLIED = "Supplied"
    SKIPPED = "Skipped"
    FAILED = "Failed"
    CANCELLED = "Cancelled"

    def __str__(self) -> str:
        return self.value

class Report(Protocol):
    """Reports back to cli.py with the current progress.
    Call it with a step name, the current status, how much is completed,
    how much total work there must be done (including completed) Leave None if
    the amount of work is not known (yet).

    Within status make the choice: Running, complete, reused, skipped or failed
    Cancelled can also be used if user canceld the run
    """

    def __call__(self, step: str, status: Status, completed: int,
                 total: int | None) -> None: ...


def init_logging(debug: bool = False,quiet: bool = False,
                 log_file: Path | None = None,
                 console: Console | None = None) -> logging.Logger:

    logger = logging.getLogger("CAT_pack")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    console = console or Console(stderr=True)

    if not quiet:
        rich_handler = RichHandler(
            console=console,
            show_time=True,
            show_path=debug,
            rich_tracebacks=debug,
            markup=True,
        )
        rich_handler.setLevel(logging.DEBUG if debug else logging.WARNING)
        logger.addHandler(rich_handler)

    if log_file is not None:
        file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s\t%(message)s"))
        logger.addHandler(file_handler)

    return logger

