import logging
from dataclasses import replace
from datetime import datetime
from decimal import Decimal
from pathlib import Path

from .options import AlignerName, BatOptions, CatOptions, ExecutionOptions, DiamondOptions, PrepareOptions
from .parsers import BinParser
from .settings import (
    BatFiles, BatSettings, CatFiles, CatSettings, ClassificationSettings, DatabaseFiles, DiamondParameters,
    DiamondSettings, ExecutionSettings, PrepareOutputs, PrepareSettings, TaxonomyFiles,
)
from .utils.check import (
    check_db_file, check_diamond, check_file, check_folder, check_integer, check_number,
    check_output_prefix, check_outputs, check_pyrodigal,
)
from .utils.errors import InputError, ValErrorCollector

log = logging.getLogger("CAT_pack")

DIAMOND_MODES = {
    "default", "faster", "fast", "mid-sensitive", "sensitive",
    "more-sensitive", "very-sensitive", "ultra-sensitive",
}


def validate_execution(options: ExecutionOptions) -> ExecutionSettings:
    return ExecutionSettings(options.threads, options.quiet, options.verbose, options.debug)


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


def validate_diamond(options: DiamondOptions) -> DiamondParameters:
    checks = ValErrorCollector()
    executable = checks.check(check_diamond, options.path_to_diamond)
    if options.mode not in DIAMOND_MODES:
        checks.add(InputError("Unknown DIAMOND mode.", hint=", ".join(DIAMOND_MODES)))
    checks.finish()
    return DiamondParameters(
        executable, options.mode, options.no_self_hits,
        options.block_size, options.index_chunks,
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
def make_prefix(args: PrepareOptions) -> PrepareOptions:
    if args.common_prefix is not None:
        return args
    return replace(args, common_prefix=f"{datetime.now():%Y-%m-%d}_CAT_pack")


def get_file_names(args: PrepareOptions) -> PrepareOutputs:
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


def validate_prepare(args: PrepareOptions) -> PrepareSettings:
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


def validate_cat_files(args: CatOptions | BatOptions, aligner: AlignerName | None) -> CatFiles | BatFiles:
    checks = ValErrorCollector()

    is_bat = isinstance(args, BatOptions)
    bins = checks.check(BinParser(args.bins, args.bin_suffix).parse) if is_bat else None
    contigs = (None if is_bat
               else checks.check(check_file, args.contigs, "Contigs file"))

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
    report = Path(f"{prefix}.{'bin' if is_bat else 'contig'}2classification.txt")

    outputs = [orf_report, report]
    proteins_gff = None
    intermediate_prefix = f"{prefix}.concatenated" if is_bat else str(prefix)
    if args.proteins is None:
        proteins = Path(f"{intermediate_prefix}.predicted_proteins.faa")
        proteins_gff = Path(f"{intermediate_prefix}.predicted_proteins.gff")
        outputs.extend((proteins, proteins_gff))
    if args.alignment is None:
        suffix = ".gz" if args.compress else ""
        alignment = Path(f"{intermediate_prefix}.alignment.{aligner or args.aligner}{suffix}")
        outputs.append(alignment)

    checks.check(check_outputs, outputs)

    checks.finish()

    report_fields = (dict(bin2contigs=bins.bin2contigs, bin_paths=bins.bin_paths)
                     if is_bat else dict(contigs=contigs))
    files_type = BatFiles if is_bat else CatFiles
    return files_type(
        proteins_fasta=proteins, proteins_gff=proteins_gff,
        alignment=alignment, fastaid2LCA=database.fastaid2LCA, branches=database.branches,
        names=database.names, nodes=database.nodes,
        orf_report=orf_report, report=report,
        **report_fields, diamond_database=database.diamond,
    )


def validate_cat(args: CatOptions | BatOptions) -> CatSettings | BatSettings:
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
            )) # TODO: Implement MMseqs2?

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

    if isinstance(args, BatOptions):
        return BatSettings(
            threads=execution.threads, quiet=execution.quiet,
            verbose=execution.verbose, debug=execution.debug,
            files=files, aligner=aligner,
            range_=classification.range_, fraction=classification.fraction,
            log_file=args.log_path, no_stars=args.no_stars,
        )
    elif isinstance(args, CatOptions):
        return CatSettings(
            threads=execution.threads, quiet=execution.quiet,
            verbose=execution.verbose, debug=execution.debug,
            files=files, aligner=aligner,
            range_=classification.range_, fraction=classification.fraction,
            log_file=args.log_path #TODO: Add no stars compatibility
        )
    else:
        raise TypeError(f"Hmmm, that type of options I don't know yet")


def get_validated_settings(args: CatOptions | BatOptions | PrepareOptions) -> CatSettings | BatSettings | PrepareSettings:
    if isinstance(args, (CatOptions, BatOptions)):
        return validate_cat(args)
    if type(args) == PrepareOptions:
        return validate_prepare(args)
    raise TypeError(f"Hmmm, that type of options class I don't know yet")


def check_orfs_match_contigs(
    contig_names: set[str],
    contig2ORFs: dict[str, list[str]],
    path: Path,
) -> None:
    overlap = len(contig_names & set(contig2ORFs))
    if overlap == 0:
        example = "contig_name_1"
        for orfs in contig2ORFs.values():
            example = orfs[0]
            break
        raise InputError(
            f"no ORFs found that can be traced back to one of the contigs "
            f"in the contigs fasta file: {example}. ORFs should be named "
            f"contig_name_#.",
            path=path,
        )

    rel_overlap = overlap / len(contig_names)
    log.info(
        f"ORFs found on {overlap:,d} / {len(contig_names):,d} contigs "
        f"({rel_overlap * 100:.2f}%)."
    )
    if rel_overlap < 0.97:
        log.warning(
            f"only {rel_overlap * 100:.2f}% contigs found with ORF predictions. This may "
            f"indicate that some contigs were missing from the protein "
            f"prediction. Please make sure that the protein prediction was "
            f"based on all contigs."
        )
