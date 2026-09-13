#!/usr/bin/env python3
import importlib

from .errors import InputError, ExternalToolError
from pathlib import Path

def check_folder(path: Path, label: str) -> Path:
    if not path.is_dir():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the folder exists.",
            path=path,
        )
    return path

def check_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the file exists."
        )


def check_db_file(folder: Path, suffix: str, label: str) -> None:
    matches = sorted(
        path for path in folder.iterdir()
        if path.is_file() and path.name.endswith(suffix)
    )

    if not matches:
        raise InputError(
            f"{label} was not found.",
            path=folder,
            hint=f"Expected a file ending in '{suffix}'.",
        )

    if len(matches) > 1:
        names = ", ".join(path.name for path in matches)

        raise InputError(
            f"Multiple files found for {label.lower()}: {names}",
            path=folder,
            hint="Use a folder containing one prepared database.",
        )

def check_pyrodigal():
    try:
        importlib.import_module("pyrodigal")
    except ImportError:
        raise ExternalToolError(
            "pyrodigal was not found",
            f"Please check whether it is installed.",
        )