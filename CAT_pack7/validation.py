from pathlib import Path
from utils.errors import InputError, CatError, ValidationError


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


def validate_args(args):
    errors = []

    def check(function, *values, **options):
        """
        This function within a function allows for a check_function
        to be called and the error to be caught and still continue to check
        other arguments and aggregate them all.

        usage: check(check_function_name, values, options)
        """
        try:
            return function(*values, **options)
        except CatError as error:
            errors.append(error)
            return None

    # check sequence inputs
    check(check_file, args.contigs, "Contigs file")

    if args.proteins is not None:
        check(check_file, args.proteins, "Protein file")

    database = check(check_folder, args.database, "Database folder")

    if database is not None:
        check(check_db_file, database,".dmnd", "DIAMOND database")
        check(check_db_file, database, "fastaid2LCAtaxid", "Prot to tax mapping")
        check(check_db_file, database, "taxids_with_multiple_offspring", "taxids_with_multiple_offspring file")

    taxonomy = check(check_folder, args.taxonomy, "Taxonomy folder")

    if taxonomy is not None:
        check(check_file, Path(taxonomy / "names.dmp"), "Taxonomy names file")
        check(check_file, Path(taxonomy / "nodes.dmp"), "Taxonomy nodes file")

    if errors:
        raise ValidationError(errors)
