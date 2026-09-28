from dataclasses import replace
from datetime import datetime
from decimal import Decimal
from pathlib import Path

from .defaults import AlignerName, CatDefaults, Defaults, DiamondDefaults, PrepareDefaults
from .settings import (
    CatFiles, CatSettings, ClassificationSettings, DatabaseFiles, DiamondParameters,
    DiamondSettings, ExecutionSettings, PrepareOutputs, PrepareSettings, TaxonomyFiles,
)
from .utils.check import (
    check_db_file, check_diamond, check_file, check_folder, check_integer, check_number,
    check_output_prefix, check_outputs, check_pyrodigal,
)
from .utils.errors import InputError, ValErrorCollector

DIAMOND_MODES = {
    "default", "faster", "fast", "mid-sensitive", "sensitive",
    "more-sensitive", "very-sensitive", "ultra-sensitive",
}


def validate_execution(defaults: Defaults) -> ExecutionSettings:
    return ExecutionSettings(defaults.threads, defaults.quiet, defaults.verbose, defaults.debug)


def validate_classification(range_, fraction) -> ClassificationSettings:
    checks = ValErrorCollector()
    range_ = checks.check(check_number, range_, "Range", 0, 100)
    fraction = checks.check(check_number, fraction, "Fraction", 0, Decimal("0.99"))
    checks.finish()
    return ClassificationSettings(range_, fraction)


def validate_taxonomy(names: Path, nodes: Path) -> TaxonomyFiles:
    checks = ValErrorCollector()
    names = checks.check(check_file, names, "Taxonomy names file")
    nodes = checks.check(check_file, nodes, "Taxonomy nodes file")
    checks.finish()
    return TaxonomyFiles(names, nodes)


def validate_aligner_name(name: str) -> AlignerName:
    normalized = name.lower()
    if normalized == "mmseqs":
        normalized = "mmseqs2"
    if normalized not in {"diamond", "mmseqs2"}:
        raise InputError("Aligner must be diamond or mmseqs2.")
    return normalized


def validate_diamond(defaults: DiamondDefaults) -> DiamondParameters:
    checks = ValErrorCollector()
    executable = checks.check(check_diamond, defaults.path_to_diamond)
    if defaults.mode not in DIAMOND_MODES:
        checks.add(InputError("Unknown DIAMOND mode.", hint=", ".join(DIAMOND_MODES)))
    checks.finish()
    return DiamondParameters(
        executable, defaults.mode, defaults.no_self_hits,
        defaults.block_size, defaults.index_chunks,
    )


def validate_database(folder: Path, *, require_diamond=False) -> DatabaseFiles:
    checks = ValErrorCollector()
    folder = checks.check(check_folder, folder, "Database folder")

    if folder is None:
        checks.finish()

    diamond = None
    if require_diamond:
        diamond = checks.check(check_db_file, folder, ".dmnd", "DIAMOND database")

    fastaid2LCA = checks.check(check_db_file, folder, "fastaid2LCAtaxid", "fastaid2LCAtaxid file")
    branches = checks.check(
        check_db_file, folder, "taxids_with_multiple_offspring",
        "taxids_with_multiple_offspring file",
    )
    taxonomy = checks.check(validate_taxonomy, folder / "names.dmp", folder / "nodes.dmp")
    checks.finish()
    return DatabaseFiles(
        fastaid2LCA=fastaid2LCA, branches=branches, diamond=diamond,
        names=taxonomy.names, nodes=taxonomy.nodes,
    )

# the replace(args) makes sure that if the run is started just before midnight
# and the get_file_names is run again after midnight the common_prefix will be
# the date of the day before
def make_prefix(args: PrepareDefaults) -> PrepareDefaults:
    if args.common_prefix is not None:
        return args
    return replace(args, common_prefix=f"{datetime.now():%Y-%m-%d}_CAT_pack")


def get_file_names(args: PrepareDefaults) -> PrepareOutputs:
    prefix = make_prefix(args).common_prefix
    folder = args.db_dir
    return PrepareOutputs(
        prefix=prefix,
        data_folder=folder,
        names=folder / "names.dmp",
        nodes=folder / "nodes.dmp",
        log_file=folder / f"{prefix}.log",
        diamond_database=folder / f"{prefix}.dmnd",
        mmseqs2_database=folder / f"{prefix}.mmseqs2",
        fastaid2LCAtaxid=folder / f"{prefix}.fastaid2LCAtaxid",
        taxids_with_multiple_offspring=folder / f"{prefix}.taxids_with_multiple_offspring",
    )


def validate_prepare(args: PrepareDefaults) -> PrepareSettings:
    args = make_prefix(args)
    checks = ValErrorCollector()
    execution = checks.check(validate_execution, args)
    db_fasta = checks.check(check_file, args.db_fasta, "Database FASTA")
    taxonomy = checks.check(validate_taxonomy, args.names, args.nodes)
    acc2tax = checks.check(check_file, args.acc2tax, "Accession-to-taxid file")
    files = get_file_names(args)

    diamond = None
    if not files.diamond_database.is_file():
        diamond = checks.check(check_diamond, args.path_to_diamond)

    checks.finish()
    return PrepareSettings(
        threads=execution.threads, quiet=execution.quiet,
        verbose=execution.verbose, debug=execution.debug,
        files=files, db_fasta=db_fasta, names=taxonomy.names,
        nodes=taxonomy.nodes, acc2tax=acc2tax, diamond=diamond, cleanup=args.cleanup,
    )


def validate_cat_files(args: CatDefaults, aligner: AlignerName | None) -> CatFiles:
    checks = ValErrorCollector()

    contigs = checks.check(check_file, args.contigs, "Contigs file")

    proteins = None
    alignment = None
    if args.proteins is not None:
        proteins = checks.check(check_file, args.proteins, "Protein file")
    elif args.alignment is None:
        checks.check(check_pyrodigal)

    if args.alignment is not None:
        alignment = checks.check(check_file, args.alignment, "Alignment file")
        if args.proteins is None:
            checks.add(InputError(
                "An existing alignment also requires its protein FASTA.",
                hint="use --proteins_fasta together with --alignment_table",
            ))

    need_diamond_db = args.alignment is None and aligner == "diamond"
    database = checks.check(validate_database, args.database, require_diamond=need_diamond_db)

    prefix = args.output_prefix
    checks.check(check_output_prefix, prefix)
    orf_report = Path(f"{prefix}.ORF2LCA.txt")
    contig_report = Path(f"{prefix}.contig2classification.txt")

    outputs = [orf_report, contig_report]
    proteins_gff = None
    if args.proteins is None:
        proteins = Path(f"{prefix}.predicted_proteins.faa")
        proteins_gff = Path(f"{prefix}.predicted_proteins.gff")
        outputs.extend((proteins, proteins_gff))
    if args.alignment is None:
        suffix = ".gz" if args.compress else ""
        alignment = Path(f"{prefix}.alignment.{aligner or args.aligner}{suffix}")
        outputs.append(alignment)

    checks.check(check_outputs, outputs)

    checks.finish()

    return CatFiles(
        contigs=contigs, proteins_fasta=proteins, proteins_gff=proteins_gff,
        alignment=alignment, fastaid2LCA=database.fastaid2LCA, branches=database.branches,
        names=database.names, nodes=database.nodes,
        orf_report=orf_report, contig_report=contig_report,
        diamond_database=database.diamond,
    )


def validate_cat(args: CatDefaults) -> CatSettings:
    checks = ValErrorCollector()
    execution = checks.check(validate_execution, args)
    classification = checks.check(validate_classification, args.range_, args.fraction)
    aligner_name = checks.check(validate_aligner_name, args.aligner)
    top = checks.check(check_integer, args.top, "top", 0, 100)
    files = checks.check(validate_cat_files, args, aligner_name)

    tmpdir = args.tmpdir if args.tmpdir is not None else args.output_prefix.parent / "tmp"
    diamond = None
    if args.alignment is None:
        if aligner_name == "diamond":
            diamond = checks.check(validate_diamond, args.diamond)
            if classification is not None and top is not None and top <= classification.range_:
                checks.add(InputError("Top must be higher than range."))
        elif aligner_name == "mmseqs2":
            checks.add(InputError(
                "MMseqs2 execution is not implemented yet.",
                hint="Use --aligner diamond or supply an existing alignment and proteins.",
            ))

    checks.finish()

    aligner = None
    if args.alignment is None:
        aligner = DiamondSettings(
            diamond=diamond.diamond, mode=diamond.mode, no_self_hits=diamond.no_self_hits,
            block_size=diamond.block_size, index_chunks=diamond.index_chunks,
            query=files.proteins_fasta, database=files.diamond_database,
            alignment=files.alignment, tmpdir=tmpdir, threads=execution.threads,
            top=top, compression=args.compress, verbose=execution.verbose
        )

    return CatSettings(
        threads=execution.threads, quiet=execution.quiet,
        verbose=execution.verbose, debug=execution.debug,
        files=files, aligner=aligner,
        range_=classification.range_, fraction=classification.fraction,
        log_file=args.log_path
    )


def get_validated_settings(args: CatDefaults | PrepareDefaults) -> CatSettings | PrepareSettings:
    if type(args) == CatDefaults:
        return validate_cat(args)
    if type(args) == PrepareDefaults:
        return validate_prepare(args)
    raise TypeError(f"Hmmm, that type of Default class I don't know yet")
