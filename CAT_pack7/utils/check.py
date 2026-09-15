#!/usr/bin/env python3
import importlib
import shutil

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

def check_file(path: Path, label: str, *, allow_empty=False) -> Path:
    if not path.is_file():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the file exists.",
            path=path,
        )

    # For files that must have data in them
    with path.open("rb") as handle:
        if not handle.read(1) and not allow_empty:
            raise InputError(f"{label} is empty", path=path,
                             hint="Remove file and try again or check if the "
                                  "correct file is supplied.")
    return path


def check_db_file(folder: Path, suffix: str, label: str, *, allow_empty=False) -> Path:
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
    return check_file(matches[0], label, allow_empty=allow_empty)

def check_pyrodigal():
    try:
        importlib.import_module("pyrodigal")
    except ImportError:
        raise ExternalToolError(
            "pyrodigal",
            "Package was not found",
            hint="Please check whether it is installed and if the correct envirnment is active",
        )

def check_diamond():
    diamon_executable = shutil.which("diamond")
    if diamon_executable is None:
        raise ExternalToolError(
            "DIAMOND",
            "was not found on PATH",
            hint="Activate the environment containing DIAMOND"
        )
    return Path(diamon_executable)
