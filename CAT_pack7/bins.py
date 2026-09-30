"""BAT!!"""
import logging
from decimal import Decimal

from . import tax
from .classification import (
    ClassificationEngine, ClassificationStatus, ORFStatus,
    import_ORFs, parse_alignment, check_orfs_match_contigs,
)
from .utils.check import check_file
from .utils.errors import InputError
from .utils.logging import Status

log = logging.getLogger("CAT_pack")


def import_bins(path, suffix):
    if path.is_dir():
        log.info(f"Importing bins from {path}.")
        paths = tuple(sorted(
            entry for entry in path.iterdir()
            if entry.is_file() and not entry.name.startswith(".")
            and entry.name.endswith(suffix) and ".concatenated." not in entry.name
        ))
        if not paths:
            raise InputError(
                f"no bins found with suffix {suffix} in bin folder. You can set the "
                "suffix with the [-s / --bin_suffix] argument.", path=path,
            )
    else:
        paths = (check_file(path, "Bin fasta"),)

    bin2contigs = {}
    contig2bin = {}
    for fasta in paths:
        contigs = bin2contigs[fasta.name] = []
        with fasta.open(encoding="utf-8") as source:
            for line in source:
                if not line.startswith(">"):
                    continue
                fields = line[1:].split()
                if not fields:
                    raise InputError("Fasta header is empty.", path=fasta)
                contig = fields[0]
                if contig in contig2bin:
                    raise InputError(
                        f"BAT has encountered {contig} twice, in {contig2bin[contig]} "
                        f"and in {fasta.name}. Fasta headers (the part "
                        "before the first space in the >line) should be unique "
                        "across bins, please remove or rename duplicates.", path=fasta,
                    )
                contig2bin[contig] = fasta.name
                contigs.append(contig)
    if not contig2bin:
        raise InputError("no contigs found in bin fasta files.", path=path)
    log.info("1 bin found!" if len(paths) == 1 else f"{len(paths):,d} bins found!")
    return bin2contigs, paths


def make_concatenated_fasta(files):
    log.info(f"Writing {files.contigs}.")
    with files.contigs.open("w", encoding="utf-8") as output:
        for path in files.bin_paths:
            with path.open(encoding="utf-8") as source:
                for line in source:
                    if line.startswith(">"):
                        output.write(f">{line[1:].split()[0]}\n")
                    else:
                        output.write(line.rstrip("\r\n") + "\n")


def bin_classification(settings, files, report):
    contig_names = {contig for contigs in files.bin2contigs.values() for contig in contigs}
    contig2orfs = import_ORFs(files.proteins_fasta)
    check_orfs_match_contigs(contig_names, contig2orfs, files.proteins_fasta)
    hits, all_hits = parse_alignment(
        files.alignment, (Decimal(100) - settings.range_) / Decimal(100),
    )
    parents, _ = tax.import_nodes(files.nodes)
    mapping = tax.import_fastaid2LCAtaxid(files.fastaid2LCA, all_hits)
    branches = tax.import_taxids_with_multiple_offspring(files.branches)
    engine = ClassificationEngine(
        taxid2parent=parents, fastaid2taxid=mapping, fraction=settings.fraction,
    )

    def lineage_text(lineage):
        lineage = list(lineage)
        if not settings.no_stars:
            lineage = tax.star_lineage(lineage, branches)
        return ";".join(lineage[::-1])

    log.info(f"BAT is flying! Files {files.bin_report} and {files.orf_report} are created.")
    classified = 0
    total = len(files.bin2contigs)
    report("Classify", Status.RUNNING, 0, total)
    with files.bin_report.open("w") as bin_out, files.orf_report.open("w") as orf_out:
        bin_out.write(
            "# bin\tclassification\treason\tlineage\t"
            f"lineage scores (f: {float(settings.fraction)})\n"
        )
        orf_out.write(
            f"# ORF\tbin\tnumber of hits (r: {settings.range_})\tlineage\ttop bit-score\n"
        )
        for index, bin_name in enumerate(sorted(files.bin2contigs), start=1):
            orfs = [orf for contig in sorted(files.bin2contigs[bin_name])
                    for orf in contig2orfs.get(contig, ())]
            result = engine.classify_group(entity_id=bin_name, orf_ids=orfs, orf2hits=hits)
            for orf in result.orf_results:
                if orf.status == ORFStatus.NO_HIT:
                    orf_out.write(f"{orf.orf_id}\t{bin_name}\tORF has no hit to database\n")
                else:
                    lineage = orf.message if orf.status == ORFStatus.NO_TAXID else lineage_text(orf.lineage)
                    orf_out.write(
                        f"{orf.orf_id}\t{bin_name}\t{orf.n_hits}\t{lineage}\t{orf.top_bitscore}\n"
                    )

            # Still trying to make this a bit more clean
            # maybe a separate writer will solve this eventually
            reasons = {
                ClassificationStatus.NO_ORFS: "no hits to database",
                ClassificationStatus.NO_HITS: "no hits to database",
                ClassificationStatus.NO_TAXIDS: "hits not found in taxonomy files",
                ClassificationStatus.NO_LINEAGE_SUPPORT: "no lineage reached minimum bit-score support",
            }
            if result.status in reasons:
                bin_out.write(f"{bin_name}\tno taxid assigned\t{reasons[result.status]}\n")
            else:
                classified += 1
                for i, assignment in enumerate(result.assignments, start=1):
                    label = "taxid assigned"
                    if len(result.assignments) > 1:
                        label += f" ({i}/{len(result.assignments)})"
                    scores = ";".join(f"{score:.2f}" for score in assignment.lineage_scores[::-1])
                    bin_out.write(
                        f"{bin_name}\t{label}\tbased on {result.based_on_n_ORFs}/{result.total_n_ORFs} ORFs\t"
                        f"{lineage_text(assignment.lineage)}\t{scores}\n"
                    )
            report("Classify", Status.RUNNING, index, total)
    log.info(
        f"BAT is done! {classified:,d}/{total:,d} bins "
        f"({classified / total * 100:.2f}%) have taxonomy assigned."
    )
    if settings.fraction < Decimal("0.5"):
        log.warning("since f is set to smaller than 0.5, one bin may have multiple classifications.")
