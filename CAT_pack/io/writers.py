"""Write the shared CAT/BAT classification report format."""
from decimal import Decimal
from importlib.resources import files
from itertools import groupby
from pathlib import Path
from typing import TextIO, IO, Literal

from CAT_pack.io.parsers import format_bitscore, FastaParser
from CAT_pack.results import Classification, ClassificationStatus, Lineage
from CAT_pack.tax import Taxonomy, ORF, convert_to_names, convert_to_official_names
from CAT_pack.utils.errors import InputError
from CAT_pack.utils.typing import EntityType


class ClassificationWriter:
    """Write headers on creation, then write the array results."""

    def __init__(
        self, classification_out: TextIO, orf_out: TextIO | None = None, *,
        entity_type: EntityType, fraction: Decimal, range_: Decimal,
        branches: set[str], no_stars: bool = False, fasta_path: Path,
        comp_bias_corr: bool = True,
    ) -> None:
        self.classification_out = classification_out
        self.orf_out = orf_out
        self.entity_type = entity_type
        self.is_bin = entity_type == "bin"
        self.fraction = fraction
        self.range_ = range_
        self.branches = branches
        self.no_stars = no_stars
        self.self_scores = self.calculate_self_bitscores(fasta_path, comp_bias_corr)
        self._lineage_taxonomy: Taxonomy | None = None
        self._lineages: dict[int, str] = {}
        self._write_headers()

    def load_self_scores(self):
        """Read only the scores of a self-hit"""
        self_scores = {}
        with files("CAT_pack").joinpath("assets/BLOSUM62").open() as f1:
            col_index = 0
            for line in f1:
                if not line.strip() or line.startswith("#"):
                    continue
                scores = line.split()
                if col_index == 0:
                    columns = scores
                else:
                    self_scores[scores[0]] = dict(zip(columns, map(int, scores[1:])))
                col_index += 1
        return self_scores

    def calculate_self_bitscores(self, proteins_fasta: Path, comp_bias_corr: bool = True) -> dict[str, Decimal]:
        self_scores = self.load_self_scores()
        background = {}
        with files("CAT_pack").joinpath("assets/BLOSUM62_background.tsv").open() as f1:
            for line in f1:
                if line.startswith("#"):
                    continue
                residue, frequency = line.split()
                background[residue] = Decimal(frequency)
        background_scores = {
            residue: sum(frequency * scores[amino_acid] for amino_acid, frequency in background.items())
            for residue, scores in self_scores.items()
        }
        log_K = Decimal("0.041").ln()
        log_2 = Decimal(2).ln()
        orf_self_bitscores = {}
        for record in FastaParser(proteins_fasta):
            sequence = record.sequence.upper().rstrip("*")
            score = Decimal(0)
            best_score = Decimal(0)
            for i, residue in enumerate(sequence):
                scores = self_scores[residue]
                correction = Decimal(0)
                if comp_bias_corr:
                    # Hauser 2016: +/-20 residues, excluding itself.
                    neighbours = sequence[max(0, i - 20):i] + sequence[i + 1:i + 21]
                    if neighbours:
                        correction = background_scores[residue] - Decimal(sum(scores[amino_acid] for amino_acid in neighbours)) / len(neighbours)
                score = max(Decimal(0), score + scores[residue] + correction)
                best_score = max(best_score, score)
            # B = (lambda * S - ln K) / ln 2. BLOSUM62, gap 11/1 table:
            # https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/src/algo/blast/core/blast_stat.c
            orf_self_bitscores[record.name] = (Decimal("0.267") * best_score - log_K) / log_2
        return orf_self_bitscores

    def _write_headers(self) -> None:
        comparison_columns = "\ttotal bitscore\tselfhit bitscore" if self.self_scores is not None else ""
        self.classification_out.write(
            f"# {self.entity_type}\tclassification\treason\tlineage\t"
            f"lineage scores (f: {float(self.fraction)}){comparison_columns}\n"
        )
        if self.orf_out:
            bin_column = "bin\t" if self.is_bin else ""
            comparison_columns = "\tselfhit bitscore" if self.self_scores is not None else ""
            self.orf_out.write(
                f"# ORF\t{bin_column}number of hits (r: {self.range_})\tlineage\ttop bit-score{comparison_columns}\n"
            )

    def write(self, result: Classification, entity: int) -> None:
        entity_id = result.entities[entity]
        orfs = result.get_orfs(entity)
        if self.orf_out is not None:
            for orf in range(orfs.start, orfs.stop):
                self.write_orf(result, orf, entity)

        status = result.get_status(entity)
        comparison = ""
        if self.self_scores is not None:
            actual = format_bitscore(int(result.weighted.total[entity]))
            self_score = sum((self.self_scores[ORF] for ORF in result.orf_names[orfs]), Decimal(0))
            comparison = f"\t{actual}\t{self_score:.1f}"
        if status != ClassificationStatus.ASSIGNED:
            empty_lineage = "\t\t" if self.self_scores is not None else ""
            self.classification_out.write(
                f"{entity_id}\tno taxid assigned\t{status.reason}{empty_lineage}{comparison}\n"
            )
            return
        if self.entity_type == EntityType.SEQUENCE:
            reason = f"based on {result.n_w_tax[orfs.start]} hits"
        else:
            reason = (
                f"based on {result.n_orfs_with_taxid[entity]}/"
                f"{orfs.stop - orfs.start} ORFs"
            )

        assignments = result.get_assignments(entity)
        weighted = result.weighted
        total = Decimal(int(weighted.total[entity]))
        score_strings: dict[int, str] = {}

        for index, assignment in enumerate(assignments, start=1):
            label = "taxid assigned"
            if len(assignments) > 1:
                label += f" ({index}/{len(assignments)})"

            start, end = weighted.lineage_offsets[assignment:assignment + 2]
            lineage = self._lineage_for(result.taxonomy, int(weighted.assignment_taxon[assignment]))
            score_parts = []
            for value in weighted.lineage_bitscores[start:end][::-1]:
                score = int(value)
                text = score_strings.get(score)
                if text is None:
                    text = f"{Decimal(score) / total:.2f}"
                    score_strings[score] = text
                score_parts.append(text)
            scores = ";".join(score_parts)
            self.classification_out.write(
                f"{entity_id}\t{label}\t{reason}\t{lineage}\t{scores}{comparison}\n"
            )

    def write_orf(self, result: Classification, orf: int, entity: int) -> None:
        prefix = result.orf_names[orf]
        if self.is_bin:
            prefix += f"\t{result.entities[entity]}"

        row = int(result.orf_row[orf])
        comparison = ""
        if self.self_scores is not None:
            comparison = f"\t{self.self_scores[result.orf_names[orf]]:.1f}"
        if row == ORF.NO_HIT:
            empty_hit = "\t\t0" if self.self_scores is not None else ""
            self.orf_out.write(f"{prefix}\tORF has no hit to database{empty_hit}{comparison}\n")
            return

        taxid = int(result.orf_lca[orf])
        if taxid == ORF.NO_TAXID:
            hits = ";".join(result.alignment.hit_names_of(row))
            lineage = f"no taxid found ({hits})"
        else:
            lineage = self._lineage_for(result.taxonomy, taxid)

        bitscore = format_bitscore(int(result.orf_top_bitscore[orf]))
        self.orf_out.write(
            f"{prefix}\t{result.orf_n_hits[orf]}\t{lineage}\t{bitscore}{comparison}\n"
        )

    def format_lineage(self, lineage: tuple[str, ...]) -> str:
        return Lineage(lineage).format(self.branches, no_stars=self.no_stars)

    def _lineage_for(self, taxonomy: Taxonomy, taxid: int) -> str:
        if taxonomy is not self._lineage_taxonomy:
            self._lineages.clear()
            self._lineage_taxonomy = taxonomy
        lineage = self._lineages.get(taxid)
        if lineage is None:
            taxids = tuple(taxonomy.taxids[position] for position in taxonomy.find_lineage(taxid))
            lineage = self.format_lineage(taxids)
            self._lineages[taxid] = lineage
        return lineage


def write_names(
    f1: TextIO,
    outf1: IO[str],
    taxid2rank: dict[str, str],
    taxid2name: dict[str, str],
    input_file: Path,
    only_official: bool = False,
    exclude_scores: bool = False,
) -> None:
    line: str | list[str]

    for line in f1:
        if line.startswith("#"):
            line = line[2:].rstrip().split("\t")

            if "lineage" not in line:
                continue

            lineage_index = line.index("lineage")

            scores_index = None
            for i in range(len(line)):
                if line[i].startswith("lineage scores"):
                    scores_index = i

            full_length = len(line)
            break
    else:
        raise InputError(f"{input_file} is not a supported classification file.", path=input_file)

    f1.seek(0)

    for line in f1:
        line = line.rstrip()

        if line.startswith("#"):
            if "lineage" not in line.lstrip("#").lstrip().split("\t"):
                outf1.write(f"{line}\n")
                continue

            if only_official:
                outf1.write(f"{line}\tsuperkingdom\tphylum\tclass\torder\t"
                            "family\tgenus\tspecies\n")
            else:
                outf1.write(f"{line}\tfull lineage names\n")

            continue

        line = line.split("\t")

        if len(line) != full_length:
            # Entry does not have a full annotation.
            # DEVNOTE: Python version 3.11 does not allow backslashes inside expression parts of f-strings
            outf1.write("{0}\n".format("\t".join(line)))
            continue

        if any([c.startswith("no taxid found") for c in line[2:4]]):
            # ORF has database hits but the accession number is not found
            # in the taxonomy files.
            outf1.write("{0}\n".format("\t".join(line)))
            continue

        if line[lineage_index] in {"", "unmapped", "unclassified"}:
            outf1.write("{0}\n".format("\t".join(line)))
            continue

        lineage = line[lineage_index].split(";")

        if scores_index and not exclude_scores:
            scores = line[scores_index].split(";")
        else:
            scores = None

        if only_official:
            names = convert_to_official_names(lineage, taxid2rank, taxid2name, scores)
        else:
            names = convert_to_names(lineage, taxid2rank, taxid2name, scores)

        outf1.write("{0}\t{1}\n".format("\t".join(line), "\t".join(names)))

def write_summary(
    f1: TextIO,
    outf1: IO[str],
    input_file: Path,
    entity_type: Literal[EntityType.CONTIG, EntityType.BIN],
    contig2length: dict[str, int],
    contigs_fasta: Path | None,
) -> None:
    is_contig = entity_type == EntityType.CONTIG
    tool = "CAT" if is_contig else "BAT"
    line: str | list[str]

    for line in f1:
        if line.startswith("#"):
            line = line.rstrip("\r\n").split("\t")

            if line[0] not in {"# contig", "# bin"}:
                continue

            if line[0] != f"# {entity_type}":
                raise InputError(
                    f"{input_file} is not a {tool} classification file.",
                    path=input_file,
                    hint=(
                        "If you want to summarise bin classifications, "
                        "simply don't supply a contigs fasta."
                        if is_contig else
                        "If you want to summarise contig classifications, "
                        "supply a contigs fasta with argument "
                        "[-c / --contigs_fasta]."
                    ),
                )

            try:
                superkingdom_index = line.index("superkingdom")
            except ValueError:
                raise InputError(
                    f"official ranks not found in header of {input_file}. "
                    f"Make sure that the {tool} classification file is named "
                    "with official ranks with "
                    "'CAT_pack add_names --only_official'.",
                    path=input_file,
                ) from None

            break
    else:
        raise InputError(
            "input file does not have a recognisable header.",
            path=input_file,
        )

    official_ranks = [
        "superkingdom", "phylum", "class", "order",
        "family", "genus", "species",
    ]

    # this replaced n_bins and inside cat there was only a n_contigs = len(contig2length)
    n_entities: dict[str, dict[str, int]] = {}
    length: dict[str, dict[str, int]] = {}
    ORFs: dict[str, dict[str, int]] = {}

    for rank in official_ranks:
        n_entities[rank] = {}
        length[rank] = {}
        ORFs[rank] = {}


    n = 0
    n_classified = 0
    total_classified_length = 0
    entity_trace = set()
    doubles = set()

    # go back to the first line
    f1.seek(0)
    for line in f1:
        line = line.rstrip("\r\n")

        if line.startswith("#"):
            continue

        n += 1
        line = line.split("\t")
        entity = line[0]

        if entity in entity_trace:
            doubles.add(entity)
        entity_trace.add(entity)

        if is_contig and entity not in contig2length:
            raise InputError(
                f"contig {entity} in CAT classification file is not found "
                "in supplied contigs fasta file. Are you sure the CAT "
                "classification file is based on the contigs fasta?",
                path=input_file,
            )

        if line[1] == "no taxid assigned":
            continue

        n_classified += 1
        ORFs_on_contig = 0

        if is_contig:
            total_classified_length += contig2length[entity]

            # NOTE that the total number of ORFs on a contig is reported,
            # not only the number of ORFs a classification is based on.
            ORFs_on_contig = int(line[2].split("/")[1].split(" ")[0])

        for i, classification in enumerate(line[superkingdom_index:]):
            classification = classification.rsplit(": ", 1)[0].rstrip("*")
            rank = official_ranks[i]

            n_entities[rank][classification] = (
                n_entities[rank].get(classification, 0) + 1
            )

            if is_contig:
                length[rank][classification] = (
                    length[rank].get(classification, 0)
                    + contig2length[entity]
                )
                ORFs[rank][classification] = (
                    ORFs[rank].get(classification, 0)
                    + ORFs_on_contig
                )

    if len(doubles) != 0:
        raise InputError(
            f"some {entity_type}s have multiple classifications. CAT_pack "
            "summarise currently does not allow for this. "
            f"{entity_type.capitalize()}s with multiple classifications: "
            f"{', '.join(list(doubles))}.",
            path=input_file,
        )

    if is_contig and n != len(contig2length):
        raise InputError(
            "the number of classified contigs is not the same as the "
            "number of contigs in contigs fasta. Are you sure the CAT "
            "classification file is based on the contigs fasta?",
            path=input_file,
        )

    if is_contig:
        n_contigs = len(contig2length)
        total_length = sum(contig2length.values())

        outf1.write(f"# total number of contigs in {contigs_fasta} is {n_contigs:,d} "
                    f"representing {total_length:,d} positions.\n")
        outf1.write(
            f"# {n_classified:,d} contigs have taxonomy assigned ({n_classified / n_contigs * 100:.2f}%) "
            f"representing {total_classified_length:,d} positions "
            f"({total_classified_length / total_length * 100:.2f}%) in {input_file}.\n"
        )

    else:
        outf1.write(
            f"# total number of bins is {n:,d}, of which {n_classified:,d} "
            f"({n_classified / n * 100:.2f}%) have taxonomy assigned.\n"
        )

    outf1.write("#\n")
    outf1.write(f"# rank\tclade\tnumber of {entity_type}s")
    if is_contig:
        outf1.write("\tnumber of ORFs\tnumber of positions")
    outf1.write("\n")

    for rank in official_ranks:
        for clade in sorted(
            n_entities[rank],
            key=lambda x: (
                length[rank][x] if is_contig else n_entities[rank][x]
            ),
            reverse=True,
        ):
            columns = [rank, clade, str(n_entities[rank][clade])]

            if is_contig:
                columns.extend([
                    str(ORFs[rank][clade]),
                    str(length[rank][clade]),
                ])

            outf1.write("\t".join(columns) + "\n")


def top_tabular_alignment(f1: TextIO, outf1: TextIO, top: int, bit_score_column: int) -> None:
    """This code manually does what DIAMOND does with the --top parameter.
    Currently, it is used for MMseqs2 alignment output.
    Sorts and filters one ORF at a time
    """
    cutoff_percent = 100 - top
    for _, lines in groupby(f1, key=lambda line: line.split("\t", 1)[0]):
        hits = [
            (int(line.rstrip().split("\t")[bit_score_column]), line)
            for line in lines
        ]
        hits.sort(key=lambda hit: hit[0], reverse=True)
        cutoff = cutoff_percent * hits[0][0]

        for bitscore, line in hits:
            if bitscore * 100 < cutoff:
                break
            outf1.write(line.rstrip("\r\n") + "\n")
