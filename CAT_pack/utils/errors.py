#!/usr/bin/env python3

# Just placeholders for now

class CatError(Exception):
    """An expected error during a CAT run."""


class InputError(CatError):
    """Invalid input files or options."""


class ExternalToolError(CatError):
    def __init__(self, tool: str, message: str):
        self.tool = tool
        super().__init__(f"{tool}: {message}")