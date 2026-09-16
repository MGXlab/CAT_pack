import logging
from pathlib import Path

from rich.console import Console
from rich.logging import RichHandler


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
