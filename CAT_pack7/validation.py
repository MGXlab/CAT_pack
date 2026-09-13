from pathlib import Path

from utils.errors import CatError, ValidationError
from utils.check import check_file, check_folder, check_db_file, check_pyrodigal


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
    else:
        check(check_pyrodigal)

    # check database folder and contents
    database = check(check_folder, args.database, "Database folder")

    if database is not None:
        check(check_db_file, database,".dmnd", "DIAMOND database")
        check(check_db_file, database, "fastaid2LCAtaxid", "Prot to tax mapping")
        check(check_db_file, database, "taxids_with_multiple_offspring", "taxids_with_multiple_offspring file")

    # check taxonomy folder and contents
    taxonomy = check(check_folder, args.taxonomy, "Taxonomy folder")

    if taxonomy is not None:
        check(check_file, Path(taxonomy / "names.dmp"), "Taxonomy names file")
        check(check_file, Path(taxonomy / "nodes.dmp"), "Taxonomy nodes file")



    if errors:
        raise ValidationError(errors)
