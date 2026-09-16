from pathlib import Path

from .tools.aligner import DiamondSettings, MMseqsSettings
from .utils.check import (
    check_file,
    check_folder,
    check_db_file,
    check_pyrodigal,
    check_diamond,
)
from .utils.errors import CatError, InputError, ValidationError


def validate_args(args):
    """Validation for the cat arguments"""
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

    file_path = {}

    file_path["contigs"] = check(check_file, args.contigs, "Contigs file")

    if args.proteins is not None:
        file_path["proteins"] = check(check_file, args.proteins, "Protein file")
    elif args.alignment is None:
        check(check_pyrodigal)

    if args.alignment is not None:
        file_path["alignment"] = check(check_file, args.alignment, "Alignment file", allow_empty=True)
        if args.proteins is None:
            errors.append(InputError("An existing alignment also requires its protein FASTA.",
                                     hint="use --proteins_fasta together with --alignment_table"))
    else:
        file_path["diamond"] = check(check_diamond)

    database = check(check_folder, args.database, "Database folder")
    if database is not None:
        if args.alignment is None:
            if args.aligner == "diamond":
                file_path["database"] = check(check_db_file, database, ".dmnd", "DIAMOND database")
        file_path["fastaid2LCA"] = check(check_db_file, database, "fastaid2LCAtaxid", "fastaid2LCAtaxid file")
        file_path["branches"] = check(check_db_file, database, "taxids_with_multiple_offspring",
                                  "taxids_with_multiple_offspring file", allow_empty=True)

    taxonomy = check(check_folder, args.taxonomy, "Taxonomy folder")
    if taxonomy is not None:
        check(check_file, taxonomy / "names.dmp", "Taxonomy names file")
        check(check_file, taxonomy / "nodes.dmp", "Taxonomy nodes file")

    prefix = args.output_prefix
    if prefix.is_dir():
        errors.append(InputError("The output prefix is a directory.", path=prefix,
                                 hint="Include a filename prefix, for example results/CAT"))
    outputs = {
        "orf_report": Path(f"{prefix}.ORF2LCA.txt"),
        "contig_report": Path(f"{prefix}.contig2classification.txt"),
    }

    if args.proteins is None:
        outputs.update(proteins_fasta=Path(f"{prefix}.predicted_proteins.faa"),
                       proteins_gff=Path(f"{prefix}.predicted_proteins.gff"))
    if args.alignment is None:
        outputs["alignment"] = Path(f"{prefix}.alignment.diamond")

    for folder in dict.fromkeys(path.parent for path in outputs.values()):
        check(check_folder, folder, "Output folder")

    path_list = []
    for path in outputs.values():
        if path.exists() or path.is_symlink():
            errors.append(InputError("Output already exists!", path=path,
                                     hint="Choose a new --output-prefix"))
        absolute_path = check(path.resolve)
        if absolute_path is None:
            continue
        if absolute_path in path_list:
            errors.append(InputError("Two outputs use the same path location", path=path))
        path_list.append(absolute_path)

    if errors:
        raise ValidationError(errors)

    file_path.update(outputs)

    return file_path


def validate_aligner_args(args, files: dict) -> DiamondSettings | MMseqsSettings | None:
    if args.alignment is not None:
        return None

    query = args.contigs
    tmpdir = args.tmpdir or args.output_prefix.parent / "tmp"
    shared = dict(
        query=query,
        database=args.database,
        alignment=args.alignment,
        tmpdir=tmpdir,
        threads=args.threads,
        compression=args.compress,
        verbose=args.verbose,
    )

    if args.aligner == "diamond":
        settings = args.diamond
        return DiamondSettings(
            diamond=files["diamond"],
            top=args.top,
            mode=settings.mode,
            no_self_hits=settings.no_self_hits,
            block_size=settings.block_size,
            index_chunks=settings.index_chunks,
            **shared,
        )

    return MMseqsSettings()
