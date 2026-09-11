#!/usr/bin/env python3
from pathlib import Path

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