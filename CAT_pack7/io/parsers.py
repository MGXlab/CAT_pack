"""input parsers"""
import gzip
import logging
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Iterator

from .. import tax
from ..config.settings import BatFiles, CatFiles
from ..utils.check import check_file
from ..utils.errors import InputError

log = logging.getLogger("CAT_pack")


@dataclass(frozen=True, slots=True)
class FastaRecord:
    name: str
    sequence: str


class FastaParser:
    """Read FASTA records and collect their contig or ORF identifiers."""

    def __init__(self, path: Path) -> None:
        self.path = path

    def __iter__(self) -> Iterator[FastaRecord]:
        # Fasta (yield) logic Adopted from:
        # https://github.com/idptools/sparrow/blob/03aa232abcbe191d96daf55ee28e1f7881979a68/sparrow/sequence_analysis/plaac/plaac.py#L318-L350
        record_name = None
        sequence_lines: list[str] = []
        with self.path.open(encoding="utf-8") as source:
            for line in source:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    next_name = self._parse_header(line)
                    if record_name is not None:
                        yield FastaRecord(record_name, "".join(sequence_lines))
                    record_name = next_name
                    sequence_lines = []
                else:
                    sequence_lines.append(line)

        if record_name is not None:
            yield FastaRecord(record_name, "".join(sequence_lines))

    def iter_headers(self) -> Iterator[str]:
        with self.path.open(encoding="utf-8") as source:
            for line in source:
                line = line.strip()
                if line.startswith(">"):
                    yield self._parse_header(line)

    def _parse_header(self, header: str) -> str:
        name = header[1:].split(" ", 1)[0]
        if not name:
            raise InputError(
                "Fasta header is empty",
                path=self.path,
                hint="Please provide a non-empty name immediately after '>'.",
            )
        return name

    def parse_contig_names(self) -> set[str]:
        """Read unique contig identifiers. Checks also for duplicate contig names"""
        log.info(f"Importing contig names from {self.path}.")
        contig_names: set[str] = set()
        for name in self.iter_headers():
            if name in contig_names:
                raise InputError(
                    "your fasta file contains duplicate headers (the "
                    "part before the first space in the >line). The "
                    f"first duplicate encountered is {name}, but there "
                    "might be more...", path=self.path,
                )
            contig_names.add(name)
        return contig_names

    def parse_ORFs(self) -> dict[str, list[str]]:
        """Group ORF identifiers by their contig_name_# prefix."""
        log.info(f"Parsing ORF file {self.path}.")
        contig2ORFs: dict[str, list[str]] = {}
        for orf_name in self.iter_headers():
            contig = orf_name.rsplit("_", 1)[0]
            if contig not in contig2ORFs:
                contig2ORFs[contig] = []
            contig2ORFs[contig].append(orf_name)
        return contig2ORFs


@dataclass(frozen=True, slots=True)
class BinInput:
    bin2contigs: dict[str, list[str]]
    bin_paths: tuple[Path, ...]


class BinParser:
    """Read bin files and reject contig names that appear more than once."""

    def __init__(self, path: Path, suffix: str) -> None:
        self.path = path
        self.suffix = suffix

    def parse(self) -> BinInput:
        paths = self._find_bin_paths()
        bin2contigs: dict[str, list[str]] = {}
        contig2bin: dict[str, str] = {}

        for path in paths:
            contigs: list[str] = []
            for contig_name in FastaParser(path).iter_headers():
                if contig_name in contig2bin:
                    previous_bin = contig2bin[contig_name]
                    raise InputError(
                        f"BAT has encountered {contig_name} twice, in {previous_bin} "
                        f"and in {path.name}. Fasta headers (the part "
                        "before the first space in the >line) should be unique "
                        "across bins, please remove or rename duplicates.",
                        path=path,
                    )
                contig2bin[contig_name] = path.name
                contigs.append(contig_name)
            bin2contigs[path.name] = contigs

        if not contig2bin:
            raise InputError("no contigs found in bin fasta files.", path=self.path)

        if len(paths) == 1:
            log.info("1 bin found!")
        else:
            log.info(f"{len(paths):,d} bins found!")
        return BinInput(bin2contigs=bin2contigs, bin_paths=paths)

    def _find_bin_paths(self) -> tuple[Path, ...]:
        if not self.path.is_dir():
            return (check_file(self.path, "Bin fasta"),)

        log.info(f"Importing bins from {self.path}.")
        paths = []
        for path in sorted(self.path.iterdir()):
            if not path.is_file():
                continue
            if path.name.startswith(".") or ".concatenated." in path.name:
                continue
            if path.name.endswith(self.suffix):
                paths.append(path)

        if not paths:
            raise InputError(
                f"no bins found with suffix {self.suffix} in bin folder. You can set the "
                "suffix with the [-s / --bin_suffix] argument.",
                path=self.path,
            )
        return tuple(paths)


@dataclass(frozen=True, slots=True)
class AlignmentInput:
    orf2hits: dict[str, list[tuple[str, Decimal]]]
    all_hits: set[str]


class AlignmentParser:
    """Keep hits within range_ percent of each ORF's highest bit score.

    Validate that input is grouped by ORF, bitscores must be descending per ORF
    """

    def __init__(self, path: Path, range_: Decimal) -> None:
        self.path = path
        self.minimum_score_fraction = (Decimal(100) - range_) / Decimal(100)

    def parse(self) -> AlignmentInput:
        log.info(f"Parsing alignment file {self.path}.")
        opener = gzip.open if self.path.suffix == ".gz" else open
        orf2hits: dict[str, list[tuple[str, Decimal]]] = {}
        all_hits: set[str] = set()
        current_orf = None
        minimum_bitscore = Decimal(0)
        below_cutoff = False
        previous_bitscore = None
        with opener(self.path, "rt", encoding="utf-8") as source:
            for line_number, line in enumerate(source, start=1):
                if not line.strip():
                    continue
                fields = line.rstrip().split("\t")
                orf_id = fields[0]
                hit_id = fields[1]
                bitscore = Decimal(fields[11])

                if orf_id != current_orf:
                    if orf_id in orf2hits:
                        raise InputError(
                            f"Alignment hits for ORF {orf_id} are not grouped at line {line_number}.",
                            path=self.path,
                            hint="Group alignment rows by ORF, with descending bit scores within each group.",
                        )
                    # The first row gives this ORF's highest bit score.
                    current_orf = orf_id
                    minimum_bitscore = self.minimum_score_fraction * bitscore
                    orf2hits[orf_id] = []
                    below_cutoff = False
                    previous_bitscore = None

                if previous_bitscore is not None and bitscore > previous_bitscore:
                    raise InputError(
                        f"Alignment bit scores for ORF {orf_id} are not descending at line {line_number}.",
                        path=self.path,
                        hint="Sort each ORF's alignment hits by descending bit score.",
                    )
                previous_bitscore = bitscore

                if below_cutoff:
                    continue
                if bitscore < minimum_bitscore:
                    # Remaining rows for this ORF have lower scores too.
                    below_cutoff = True
                    continue

                orf2hits[orf_id].append((hit_id, bitscore))
                all_hits.add(hit_id)

        return AlignmentInput(orf2hits=orf2hits, all_hits=all_hits)


@dataclass(frozen=True, slots=True)
class ClassificationInputs:
    entity_type: str
    entity2ORFs: dict[str, list[str]]
    contig_names: set[str]
    contig2ORFs: dict[str, list[str]]
    alignment: AlignmentInput
    taxid2parent: dict[str, str]
    fastaid2taxid: dict[str, str]
    branches: set[str]


class ClassificationParser:
    """Load ORF groups, accepted alignment hits, and their taxonomy data."""

    def __init__(self, files: CatFiles | BatFiles, range_: Decimal) -> None:
        self.files = files
        self.range_ = range_

    def parse(self) -> ClassificationInputs:
        contig2ORFs = FastaParser(self.files.proteins_fasta).parse_ORFs()
        if isinstance(self.files, BatFiles):
            entity_type = "bin"
            contig_names = set()
            for contigs in self.files.bin2contigs.values():
                contig_names.update(contigs)
            entity2ORFs = self._group_bin_ORFs(self.files.bin2contigs, contig2ORFs)
        else:
            entity_type = "contig"
            contig_names = FastaParser(self.files.contigs).parse_contig_names()
            entity2ORFs = {}
            for contig in contig_names:
                entity2ORFs[contig] = contig2ORFs.get(contig, [])

        alignment = AlignmentParser(self.files.alignment, self.range_).parse()
        taxid2parent, _ = tax.import_nodes(self.files.nodes)
        fastaid2taxid = tax.import_fastaid2LCAtaxid(self.files.fastaid2LCA, alignment.all_hits)
        branches = tax.import_taxids_with_multiple_offspring(self.files.branches)

        return ClassificationInputs(
            entity_type=entity_type,
            entity2ORFs=entity2ORFs,
            contig_names=contig_names,
            contig2ORFs=contig2ORFs,
            alignment=alignment,
            taxid2parent=taxid2parent,
            fastaid2taxid=fastaid2taxid,
            branches=branches,
        )

    def _group_bin_ORFs(
        self,
        bin2contigs: dict[str, list[str]],
        contig2ORFs: dict[str, list[str]],
    ) -> dict[str, list[str]]:
        bin2ORFs: dict[str, list[str]] = {}
        for bin_name, contigs in bin2contigs.items():
            orfs: list[str] = []
            for contig in sorted(contigs):
                orfs.extend(contig2ORFs.get(contig, []))
            bin2ORFs[bin_name] = orfs
        return bin2ORFs
