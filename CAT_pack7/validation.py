from pathlib import Path
from utils.errors import InputError, CatError, ValidationError


def check_folder(path: Path, label: str) -> None:
    if not path.is_dir():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the folder exists.",
            path=path,
        )

def check_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise InputError(
            f"{label} was not found.",
            hint="Double check if the file exists."
        )



def validate_args(args):
    errors = []

    def check(function, *values, **options):
        """Run a check, remember expected errors, and continue."""
        try:
            return function(*values, **options)
        except CatError as error:
            errors.append(error)
            return None

    # check sequence inputs
    check(check_file, args.contigs, "Contigs file")

    if args.proteins is not None:
        check(check_file, args.proteins, "Protein file")

    check(check_folder, args.database, "Database folder")
    check(check_folder, args.taxonomy, "Taxonomy folder")

    if errors:
        raise ValidationError(errors)
