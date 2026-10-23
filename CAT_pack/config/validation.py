import logging
import shutil
import sys
from dataclasses import replace, fields
from datetime import datetime
from decimal import Decimal
from pathlib import Path

from .options import AlignerName, BatOptions, CatOptions, ExecutionOptions, DiamondOptions, PrepareOptions, \
    MemoryOptions, DownloadOptions, RatOptions, AnnotationOptions
from .settings import (
    BatFiles, BatSettings, CatFiles, CatSettings, ClassificationSettings, DatabaseFiles, DiamondParameters,
    DiamondPrepareSettings, DiamondSettings, ExecutionSettings, MMseqsPrepareSettings,
    PrepareOutputs, PrepareSettings, TaxonomyFiles, MemorySettings,
    DownloadFile, DownloadOutputs, DownloadSettings, NRUrls, GTDBUrls, GTDBDownloadOutputs, RatFiles, RatSettings,
)
from ..io.parsers import BinParser
from ..utils.check import (
    check_db_file, check_diamond, check_fasta_file, check_folder, check_integer, check_number,
    check_output_prefix, check_outputs, check_pyrodigal, check_file_not_empty, check_file,
)
from ..utils.errors import ExternalToolError, InputError, ValErrorCollector
from ..utils.memory import GIB

log = logging.getLogger("CAT_pack")

DIAMOND_MODES = {
    "default", "faster", "fast", "mid-sensitive", "sensitive",
    "more-sensitive", "very-sensitive", "ultra-sensitive",
}


def validate_execution(options: ExecutionOptions) -> ExecutionSettings:
    memory = validate_memory(options.memory) if isinstance(options, (CatOptions, BatOptions, RatOptions)) else MemorySettings()
    if memory.available_memory_bytes is not None:
        log.info(f"[red]Memory mode[/red]: custom ({memory.available_memory_bytes / GIB:g} GiB)")
    elif memory.low_memory:
        log.info("Memory mode: low")
    else:
        log.info("[red]Memory mode[/red]: automatic")
    return ExecutionSettings(options.threads, options.quiet, options.verbose, options.debug, memory=memory)


def validate_classification(range_, fraction) -> ClassificationSettings:
    checks = ValErrorCollector()
    log.debug("Checking range parameter")
    range_ = checks.check(check_number, range_, "Range", 0, 100)
    log.info(f"Checked Range parameter: {range_}")
    log.debug("Checking fraction parameter")
    fraction = checks.check(check_number, fraction, "Fraction", 0, Decimal("0.99"))
    log.info(f"Checked Fraction parameter: {fraction}")
    checks.finish()
    return ClassificationSettings(range_, fraction)


def validate_memory(options: MemoryOptions) -> MemorySettings:
    if sum((options.low_memory, options.high_memory, options.available_memory is not None)) > 1:
        raise InputError("Woah, too many memories",
                         hint="Choose only one of --low-memory, --high-memory, --available-memory.")
    budget = None
    if options.available_memory is not None:
        budget = check_integer(options.available_memory, "Available memory (GiB)", 1, sys.maxsize) * GIB
    return MemorySettings(low_memory=options.low_memory, available_memory_bytes=budget)


def validate_taxonomy(names: Path, nodes: Path) -> TaxonomyFiles:
    checks = ValErrorCollector()
    names = checks.check(check_file_not_empty, names, "Taxonomy names file")
    nodes = checks.check(check_file_not_empty, nodes, "Taxonomy nodes file")
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


def get_download_file_names(args: DownloadOptions) -> DownloadOutputs:
    prefix = f"{datetime.now():%Y-%m-%d}"
    folder = args.out_dir
    db_type = args.db_type.lower().strip()
    if db_type not in {"nr", "gtdb"}:
        raise InputError("Incorrect database selected. Please select nr or gtdb.")
    if db_type == "gtdb":
        urls = {url.rsplit("/", 1)[-1]: url for url in GTDBUrls().urls}
        paths = {name: folder / f"{prefix}.{name}" for name in urls}
        return GTDBDownloadOutputs(
            prefix=prefix, data_folder=folder,
            log_file=folder / f"{prefix}.download.log",
            db_fasta=folder / f"{prefix}.gtdb_seqs.fa.gz",
            acc2tax=folder / f"{prefix}.prot.accession2taxid.txt.gz",
            taxonomy_archive=None,
            names=folder / f"{prefix}.names.dmp",
            nodes=folder / f"{prefix}.nodes.dmp",
            downloads=tuple(DownloadFile(url, paths[name]) for name, url in urls.items()),
            protein_archive=paths["gtdb_proteins_aa_reps.tar.gz"],
            taxonomy_tables=(paths["ar53_taxonomy.tsv.gz"], paths["bac120_taxonomy.tsv.gz"]),
            trees=(paths["ar53.tree"], paths["bac120.tree"]),
            tree=folder / f"{prefix}.gtdb.tree",
            version=paths["VERSION.txt"], checksums=paths["MD5SUM.txt"],
        )
    urls = {url.rsplit("/", 1)[-1]: url for url in NRUrls().urls}
    paths = {name: folder / f"{prefix}.{name}" for name in urls}
    downloads = []
    for name, url in urls.items():
        downloads.append(DownloadFile(url, paths[name], paths.get(f"{name}.md5")))
    return DownloadOutputs(
        prefix=prefix, data_folder=folder,
        log_file=folder / f"{prefix}.download.log",
        db_fasta=paths["nr.gz"],
        acc2tax=paths["prot.accession2taxid.FULL.gz"],
        taxonomy_archive=paths["taxdump.tar.gz"],
        names=folder / f"{prefix}.names.dmp",
        nodes=folder / f"{prefix}.nodes.dmp",
        downloads=tuple(downloads),
    )


def validate_download(args: DownloadOptions, files: DownloadOutputs | None = None) -> DownloadSettings:
    match args.db_type.lower().strip():
        case "nr":
            db_type = "nr"
        case "gtdb":
            db_type = "gtdb"
        case _:
            raise InputError("Incorrect database selected. Please select nr or gtdb.")
    files = files or get_download_file_names(args)
    checks = ValErrorCollector()
    if files.data_folder.exists() and not files.data_folder.is_dir():
        checks.add(InputError("Download output location must be a directory.", path=files.data_folder))
    outputs = {item.path for item in files.downloads} | {
        files.names, files.nodes, files.db_fasta, files.acc2tax,
    }
    if not args.no_log:
        outputs.add(files.log_file)
    if isinstance(files, GTDBDownloadOutputs):
        outputs.add(files.tree)
    for path in outputs:
        checks.check(check_file, path, path.name, exist=False)
    checks.finish()

    return DownloadSettings(db_type=db_type, files=files, cleanup=args.cleanup,
                        quiet=args.quiet, no_log=args.no_log, debug=args.debug)


def validate_prepare(args: PrepareOptions, files: PrepareOutputs | None = None) -> PrepareSettings:
    args = make_prefix(args)
    checks = ValErrorCollector()
    execution = checks.check(validate_execution, args)
    db_fasta = checks.check(check_fasta_file, args.db_fasta, "Database FASTA")
    taxonomy = checks.check(validate_taxonomy, args.names, args.nodes)
    acc2tax = checks.check(check_file_not_empty, args.acc2tax, "Accession-to-taxid file")
    files = files or get_file_names(args)
    threads = checks.check(check_integer, args.threads, "Threads", 1, sys.maxsize)
    if not files.prefix or files.prefix in {".", ".."} or any(c in files.prefix for c in '/\\:'):
        checks.add(InputError("The prefix contains directory separators",
                              hint="The prefix must be a filename without directory separators."))
    if files.data_folder.exists() and not files.data_folder.is_dir():
        checks.add(InputError("Database output location must be a directory.", path=files.data_folder))
    for output in (files.names, files.nodes, files.log_file, files.diamond_database,
                   files.mmseqs2_database, files.fastaid2LCAtaxid,
                   files.taxids_with_multiple_offspring):
        if output.exists() and not output.is_file():
            checks.add(InputError("Expected a file at the output location.", path=output))

    diamond = None
    if not files.diamond_database.is_file():
        diamond = checks.check(check_diamond, args.path_to_diamond)

    # Keeps mmseq2 use optional
    mmseqs = None
    if args.build_mmseqs2 and not files.mmseqs2_database.is_file():
        command = str(args.path_to_mmseqs.resolve()) if args.path_to_mmseqs else "mmseqs"
        found = shutil.which(command)
        if found is None:
            checks.add(ExternalToolError("MMseqs2", "was not found", hint="Supply --path-to-mmseqs."))
        else:
            mmseqs = Path(found)

    checks.finish()
    return PrepareSettings(
        threads=threads, quiet=execution.quiet,
        verbose=execution.verbose, debug=execution.debug,
        files=files, db_fasta=db_fasta, names=taxonomy.names,
        nodes=taxonomy.nodes, acc2tax=acc2tax,
        diamond=(DiamondPrepareSettings(diamond, db_fasta, files.diamond_database,
                                       threads, execution.verbose) if diamond is not None else None),
        mmseqs=(MMseqsPrepareSettings(mmseqs, db_fasta, files.mmseqs2_database,
                                     threads, execution.verbose) if mmseqs is not None else None),
    )


def validate_cat_files(args: CatOptions | BatOptions, aligner: AlignerName | None) -> CatFiles | BatFiles:
    checks = ValErrorCollector()

    is_bat = isinstance(args, BatOptions)
    bins = None
    if is_bat:
        parser = BinParser(args.bins, args.bin_suffix)
        bin_paths = checks.check(parser.find_bin_paths)
        if bin_paths is not None:
            valid_paths = [checks.check(check_fasta_file, path, "Bin fasta") for path in bin_paths]
            if all(path is not None for path in valid_paths):
                bins = checks.check(parser.parse, bin_paths)
    contigs = (None if is_bat
               else checks.check(check_fasta_file, args.contigs, "Contigs file"))

    proteins = None
    alignment = None
    if args.proteins is not None:
        proteins = checks.check(check_fasta_file, args.proteins, "Protein file")
    elif args.alignment is None:
        checks.check(check_pyrodigal)

    if args.alignment is not None:
        alignment = checks.check(check_file_not_empty, args.alignment, "Alignment file")
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
    log.info(f"Validating execution arguments")
    execution = checks.check(validate_execution, args)
    log.info(f"Validating classification arguments")
    classification = checks.check(validate_classification, args.range_, args.fraction)
    log.info(f"Validating alignment arguments")
    aligner_name = checks.check(validate_aligner_name, args.aligner)
    top = checks.check(check_integer, args.top, "top", 0, 100)
    log.info(f"Validating file arguments")
    files = checks.check(validate_cat_files, args, aligner_name)

    tmpdir = args.tmpdir if args.tmpdir is not None else args.output_prefix.parent / "tmp"
    diamond = None
    if args.alignment is None:
        if aligner_name == "diamond":
            log.debug(f"Using diamond")
            log.info(f"Validating diamond arguments")
            diamond = checks.check(validate_diamond, args.diamond)
            if classification is not None and top is not None and top <= classification.range_:
                checks.add(InputError("Top must be higher than range."))
        elif aligner_name == "mmseqs2":
            log.debug(f"Using mmseqs2")
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
            memory=execution.memory,
            range_=classification.range_, fraction=classification.fraction,
            log_file=args.log_path, no_stars=args.no_stars,
        )
    elif isinstance(args, CatOptions):
        return CatSettings(
            threads=execution.threads, quiet=execution.quiet,
            verbose=execution.verbose, debug=execution.debug,
            files=files, aligner=aligner,
            memory=execution.memory,
            range_=classification.range_, fraction=classification.fraction,
            log_file=args.log_path #TODO: Add no stars compatibility
        )
    else:
        raise TypeError(f"Hmmm, that type of options I don't know yet")


def validate_rat_files(args: RatOptions) -> RatFiles:
    checks = ValErrorCollector()
    prefix = str(args.output_prefix)
    output_files = []
    checks.check(check_fasta_file, args.contigs_fasta, "Contigs FASTA")
    nodes = checks.check(check_file_not_empty, args.database / "nodes.dmp", "Taxonomy nodes file")
    checks.check(check_output_prefix, args.output_prefix)
    checks.check(check_output_prefix, args.log_path)
    checks.check(check_integer, args.threads, "Threads", 1, sys.maxsize)
    checks.check(check_integer, args.mapping_quality, "Mapping quality", 0, 255)
    if args.sort_memory_mib is not None:
        checks.check(check_integer, args.sort_memory_mib, "Sort memory (MiB)", 1, sys.maxsize)
    checks.check(validate_classification, args.range_, args.fraction)
    checks.check(validate_memory, args.memory)


    # Check the optional arguments and argument combinations
    read_files = []
    if args.read_file is not None:
        if len(args.read_file) > 2:
            checks.add(InputError("Too many read files", hint="supply one or two read files."))
        for read_file in args.read_file:
            read_files.append(checks.check(check_file_not_empty, read_file, "Read file"))

    mapping_sam = None
    mapping_bam = None
    bam_files = ()
    if args.bam_file is not None:
        if len(args.bam_file) > 2:
            checks.add(InputError("Too many bam files", hint="supply one or two bam files. One is preferred"))
        for bam_file in args.bam_file:
            checks.check(check_file_not_empty, bam_file, "Bam file")
        bam_files = args.bam_file
    elif args.read_file:
        mapping_sam = Path(f"{prefix}.{args.contigs_fasta.name}.{args.read_file[0].name}.bwamem")
        mapping_bam = Path(f"{mapping_sam}.bam")
        bam_files = (Path(f"{mapping_sam}.sorted"),)
        output_files += [mapping_sam, mapping_bam, *bam_files]

    if not args.read_file and not args.bam_file:
        checks.add(InputError("No read files or bam files supplied.", hint="supply read files (--read_file) or bam files (--bam_file)."))

    if args.contigs2classification is not None:
        c2c = checks.check(check_file_not_empty, args.contigs2classification, "Contigs classification")
    else:
        c2c = Path(f"{prefix}.CAT.contig2classification.txt")
    if args.bin2classification is not None:
        checks.check(check_file_not_empty, args.bin2classification, "Bin classification")
    if args.unmapped2classification is not None:
        checks.check(check_file_not_empty, args.unmapped2classification, "Unmapped classification")
    if args.proteins is not None:
        checks.check(check_fasta_file, args.proteins, "Protein FASTA")
    if args.alignment is not None:
        checks.check(check_file_not_empty, args.alignment, "Alignment table")
    if args.alignment_unmapped is not None:
        checks.check(check_file, args.alignment_unmapped, "Direct alignment")

    bins = None
    if "m" in args.mode.value:
        if args.bin_fasta_or_folder is None:
            checks.add(InputError("No bins supplied", hint="Add the require bins (--bin_fasta or --bin-folder)."))
        else:
            parser = BinParser(args.bin_fasta_or_folder, args.bin_suffix)
            paths = checks.check(parser.find_bin_paths)
            if paths is not None:
                for path in paths:
                    checks.check(check_fasta_file, path, "Bin FASTA")
                bins = checks.check(parser.parse, paths)
    direct = "r" in args.mode.value
    files = RatFiles(
        contigs_fasta=args.contigs_fasta,
        nodes=nodes,
        read_files=tuple(read_files),
        bam_files=bam_files,
        mapping_sam=mapping_sam,
        mapping_bam=mapping_bam,
        bin_paths=bins.bin_paths if bins else (),
        bin2contigs=bins.bin2contigs if bins else {},
        contigs2classification=c2c if "c" in args.mode.value else None,
        bin2classification=(args.bin2classification or Path(f"{prefix}.BAT.bin2classification.txt")) if "m" in args.mode.value else None,
        unmapped2classification=(args.unmapped2classification or Path(f"{prefix}.RAT.unmapped2classification.txt")) if direct else None,
        unclassified_fasta=Path(f"{prefix}.RAT.unclassified.fna") if direct else None,
        alignment_unmapped=(args.alignment_unmapped or Path(f"{prefix}.RAT.unclassified.alignment.diamond")) if direct else None,
        contig_abundance=Path(f"{prefix}.RAT.contig.abundance.txt"),
        complete_abundance=Path(f"{prefix}.RAT.complete.abundance.txt"),
        bin_abundance=Path(f"{prefix}.RAT.bin.reads.txt") if bins else None,
        read_classification=Path(f"{prefix}.RAT.read2classification.txt") if args.read2classification else None,
    )
    output_files += [files.contig_abundance, files.complete_abundance]
    if files.bin_abundance is not None: output_files.append(files.bin_abundance)
    if files.read_classification is not None: output_files.append(files.read_classification)
    if "r" in args.mode.value and not args.unmapped2classification:
        output_files.append(files.unmapped2classification)
        if not args.alignment_unmapped:
            output_files.append(files.unclassified_fasta)
            output_files.append(files.alignment_unmapped)

    checks.finish()
    return replace(files, outputs=tuple(output_files))



def make_cat_options(args: RatOptions) -> CatOptions:
    shared = {
        field.name: getattr(args, field.name)
        for field in fields(AnnotationOptions)
    }

    shared.update(
        output_prefix=Path(f"{args.output_prefix}.CAT"),
        log_file=args.log_path,
    )

    return CatOptions(
        contigs=args.contigs_fasta,
        **shared,
    )

def make_bat_options(args: RatOptions) -> BatOptions:
    shared = {
        field.name: getattr(args, field.name)
        for field in fields(AnnotationOptions)
    }

    shared.update(
        output_prefix=Path(f"{args.output_prefix}.BAT"),
        log_file=args.log_path
    )

    return BatOptions(bins=args.bin_fasta_or_folder, bin_suffix=args.bin_suffix,
                      no_stars=args.no_stars, **shared)


def validate_rat(args: RatOptions) -> RatSettings:
    checks = ValErrorCollector()
    files = checks.check(validate_rat_files, args)
    execution = checks.check(validate_execution, args)
    classification = checks.check(validate_classification, args.range_, args.fraction)
    bwa = checks.check(shutil.which, str(args.path_to_bwa))
    samtools = checks.check(shutil.which, str(args.path_to_samtools))

    if bwa is None:
        checks.add(ExternalToolError("BWA","executable not found."))
    if samtools is None:
        checks.add(ExternalToolError("SAMtools","executable not found."))

    cat_settings = None
    if "c" in args.mode.value and args.contigs2classification is None:
        cat_settings = checks.check(validate_cat, make_cat_options(args))

    bat_settings = None
    if "m" in args.mode.value and args.bin2classification is None:
        bat_settings = checks.check(validate_cat, make_bat_options(args))

    aligner = None
    if "r" in args.mode.value and args.alignment_unmapped is None and args.unmapped2classification is None:
        diamond = checks.check(validate_diamond, args.diamond)
        database = checks.check(
            check_db_file, args.database, ".dmnd", "DIAMOND database",
        )
        top = checks.check(check_integer, args.top, "Top", 0, 100)

        if top <= classification.range_:
            checks.add(InputError("Top must be higher than range."))

        #checks.finish before the aligner arguments are loaded t
        checks.finish()
        aligner = DiamondSettings(
            diamond=diamond.diamond, mode=diamond.mode,
            no_self_hits=diamond.no_self_hits, block_size=diamond.block_size,
            index_chunks=diamond.index_chunks, query=files.unclassified_fasta,
            database=database, alignment=files.alignment_unmapped,
            tmpdir=(args.tmpdir or args.output_prefix.parent)
                   / f"{args.output_prefix.name}.RAT.tmp",
            threads=execution.threads, top=top, compression=args.compress,
            verbose=execution.verbose, blast_flavour="blastx", tiny_out=True,
        )
    else:
        checks.finish()


    return RatSettings(
        files=files, log_file=args.log_path,
        bwa=Path(bwa), samtools=Path(samtools),
        sort_memory_mib=args.sort_memory_mib,
        threads=execution.threads, quiet=execution.quiet,
        verbose=execution.verbose, debug=execution.debug,
        memory=execution.memory, cat=cat_settings, bat=bat_settings, mapping_quality=args.mapping_quality, aligner=aligner
    )

def get_validated_settings(args: CatOptions | BatOptions | RatOptions | PrepareOptions) -> CatSettings | BatSettings | RatSettings | PrepareSettings:
    if isinstance(args, RatOptions):
        return validate_rat(args)
    if isinstance(args, (CatOptions, BatOptions)):
        log.debug(f"Starting validation on CatOptions or BatOptions")
        return validate_cat(args)
    if type(args) == PrepareOptions:
        log.debug(f"Starting validation on PrepareOptions")
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
