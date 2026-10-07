#!/usr/bin/env python3
import importlib
import os
import shutil
from decimal import Decimal, InvalidOperation
from pathlib import Path

from .errors import InputError, ExternalToolError, ValErrorCollector


def check_folder(path: Path, label: str) -> Path:
    if not path.is_dir():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the folder exists.",
            path=path,
        )
    return path


def check_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the file exists.",
            path=path,
        )
    return path

def check_fasta_file(path: Path, label: str) -> Path:
    path = check_file(path, label)
    # check if fasta file is empty
    os.lockf(path)


def check_db_file(folder: Path, suffix: str, label: str) -> Path:
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
    return matches[0]


def check_pyrodigal():
    try:
        importlib.import_module("pyrodigal")
    except ImportError:
        raise ExternalToolError(
            "pyrodigal",
            "Package was not found",
            hint="Please check whether it is installed and if the correct environment is active",
        )


def check_diamond(path: Path | None = None) -> Path:
    command = str(path.resolve()) if path is not None else "diamond"
    found = shutil.which(command)
    if found is None:
        raise ExternalToolError(
            "DIAMOND",
            "was not found on PATH",
            hint="Activate the environment containing DIAMOND",
        )
    return Path(found)


def check_number(value, label: str, minimum, maximum) -> Decimal:
    try:
        number = Decimal(str(value))
    except (InvalidOperation, ValueError):
        number = None
    if number is not None and number.is_finite() and minimum <= number <= maximum:
        return number
    raise InputError(f"{label} must be a number between {minimum} and {maximum}.")


def check_integer(value, label: str, minimum, maximum) -> int:
    number = check_number(value, label, minimum, maximum)
    if number != number.to_integral_value():
        raise InputError(f"{label} must be an integer between {minimum} and {maximum}.")
    return int(number)


def check_output_prefix(prefix: Path) -> Path:
    if prefix.is_dir():
        raise InputError(
            "prefix for output files is a directory.",
            path=prefix,
            hint="Include a filename prefix, for example results/CAT",
        )
    if not prefix.parent.is_dir():
        raise InputError(
            f"cannot find output directory {prefix.parent} "
                    f"to which output files should be written.",
            path=prefix.parent,
        )
    return prefix


def check_outputs(paths: list[Path]) -> None:
    checks = ValErrorCollector()
    for path in paths:
        if path.is_file():
            checks.add(InputError(
                "Output already exists!",
                path=path,
                hint="Choose a different output location",
            ))
    checks.finish()
