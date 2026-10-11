"""input parsers"""
import bz2
import gzip
import logging
from array import array
from dataclasses import dataclass
from decimal import Decimal
from itertools import chain
from pathlib import Path
from typing import Iterator

import numpy as np

from .reference_mapping import read_fastaid2LCAtaxid
from .taxonomy import read_branch_taxids
from .. import tax
from ..config.settings import BatFiles, CatFiles, ParserSettings
from ..utils.errors import InputError

log = logging.getLogger("CAT_pack")
parser_settings = ParserSettings()


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
        opener = gzip.open if self.path.suffix == ".gz" else bz2.open if self.path.suffix == ".bz2" else open
        with opener(self.path, "rt", encoding="utf-8") as source:
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
        opener = gzip.open if self.path.suffix == ".gz" else bz2.open if self.path.suffix == ".bz2" else open
        with opener(self.path, "rt", encoding="utf-8") as source:
            for line in source:
                line = line.lstrip()
                if line.startswith(">"):
                    yield self._parse_header(line.rstrip())

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

    def parse(self, paths=None) -> BinInput:
        paths = self.find_bin_paths() if paths is None else paths
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

    def find_bin_paths(self) -> tuple[Path, ...]:
        if not self.path.is_dir():
            return (self.path,)

        log.info(f"Importing bins from {self.path}.")
        paths = []
        for path in sorted(self.path.iterdir()):
            if not path.is_file():
                continue
            if path.name.startswith(".") or ".concatenated." in path.name:
                continue
            if path.name.endswith((self.suffix, f"{self.suffix}.gz", f"{self.suffix}.bz2")):
                paths.append(path)

        if not paths:
            raise InputError(
                f"no bins found with suffix {self.suffix} in bin folder. You can set the "
                "suffix with the [-s / --bin_suffix] argument.",
                path=self.path,
            )
        return tuple(paths)


def parse_bitscore(text: str) -> int:
    whole_number, _, decimal = text.strip().partition(".")
    return int(whole_number) * parser_settings.bit_score_multiplier + int(decimal or "0")


def format_bitscore(bitscore: int) -> str:
    whole_number = bitscore // parser_settings.bit_score_multiplier
    decimal = bitscore % parser_settings.bit_score_multiplier
    if decimal == 0:
        return str(whole_number)
    return f"{whole_number}.{decimal}"

@dataclass(frozen=True, slots=True)
class AlignmentInput:
    orf_names: list[str]
    orf_index: dict[str, int]
    orf_offsets: np.ndarray
    hit_ids: np.ndarray
    bitscores: np.ndarray
    hit_names: list[str]

    def hits_of(self, orf: int) -> slice:
        return slice(self.orf_offsets[orf], self.orf_offsets[orf + 1])

    def hit_names_of(self, orf: int) -> list[str]:
        return [self.hit_names[hit] for hit in self.hit_ids[self.hits_of(orf)].tolist()]


class AlignmentParser:
    """Keep hits within range_ percent of each ORF's highest bit score.

    Validate that input is grouped by ORF, bitscores must be descending per ORF
    """

    def __init__(self, path: Path, range_: Decimal, cat_only: bool | None = None) -> None:
        self.path = path
        self.cutoff_numerator, self.cutoff_denominator = ((Decimal(100) - range_)/Decimal(100)).as_integer_ratio()
        self.cat_only = cat_only

    def parse(self) -> AlignmentInput:
        log.info(f"Parsing alignment file {self.path}.")
        if self.path.suffix == ".gz":
            source_file = gzip.open(self.path, "rt", encoding="utf-8")
        elif self.path.suffix == ".bz2":
            source_file = bz2.open(self.path, "rt", encoding="utf-8")
        else:
            source_file = open(self.path, "rt", encoding="utf-8", buffering=1024 * 1024)
        orf_names: list[str] = []
        orf_index: dict[str, int] = {}
        orf_offsets = array("q")
        hit_ids = array("i")
        bitscores = array("q")
        hit2index: dict[str, int] = {}
        with source_file as source:
            current_orf = None
            cutoff = 0
            below_cutoff = False
            previous_bitscore = None
            bitscore_values: dict[str, int] = {}
            lines = enumerate(source, start=1)
            first_row = next((row for row in lines if row[1].strip()), None)
            rows = chain((first_row,), lines) if first_row is not None else ()
            if self.cat_only is None:
                bitscore_index = 11
            else:
                bitscore_index = 2 if self.cat_only else 11

            for line_number, line in rows:
                if below_cutoff:
                    first_tab = line.find("\t")
                    if first_tab != -1 and line[:first_tab] == current_orf:
                        continue
                if line == "\n" or line == "\r\n":
                    continue
                fields = line.split("\t", bitscore_index + 1)
                orf_id = fields[0]
                hit_id = fields[1]
                bitscore_string = fields[bitscore_index].rstrip("\r\n")
                bitscore = bitscore_values.get(bitscore_string)
                if bitscore is None:
                    bitscore = parse_bitscore(bitscore_string)
                    bitscore_values[bitscore_string] = bitscore

                if orf_id != current_orf:
                    if orf_id in orf_index:
                        raise InputError(
                            f"Alignment hits for ORF {orf_id} are not grouped at line {line_number}.",
                            path=self.path,
                            hint="Group alignment rows by ORF, with descending bit scores within each group.")
                    current_orf = orf_id
                    cutoff = bitscore * self.cutoff_numerator
                    orf_index[orf_id] = len(orf_names)
                    orf_names.append(orf_id)
                    orf_offsets.append(len(hit_ids))
                    below_cutoff = False
                    previous_bitscore = None

                if previous_bitscore is not None and bitscore > previous_bitscore:
                    raise InputError(
                        f"Alignment bit scores for ORF {orf_id} are not descending at line {line_number}.",
                        path=self.path,
                        hint="Sort each ORF's alignment hits by descending bit score.")
                previous_bitscore = bitscore

                if bitscore * self.cutoff_denominator < cutoff:
                    below_cutoff = True
                    continue

                hit_ids.append(hit2index.setdefault(hit_id, len(hit2index)))
                bitscores.append(bitscore)

        orf_offsets.append(len(hit_ids))
        return AlignmentInput(
            orf_names=orf_names,
            orf_index=orf_index,
            orf_offsets=np.frombuffer(orf_offsets, dtype=np.int64),
            hit_ids=np.frombuffer(hit_ids, dtype=np.int32),
            bitscores=np.frombuffer(bitscores, dtype=np.int64),
            hit_names=list(hit2index),
        )


@dataclass(frozen=True, slots=True)
class ClassificationInputs:
    entity_type: str
    entity2ORFs: dict[str, list[str]]
    contig_names: set[str]
    contig2ORFs: dict[str, list[str]]
    alignment: AlignmentInput
    taxonomy: tax.Taxonomy
    fastaid2taxid: np.ndarray  # tax.NO_TAXID for absent or unknown taxids
    branches: set[str]


class ClassificationParser:
    """Load ORF groups, accepted alignment hits, and their taxonomy data."""

    def __init__(self, files: CatFiles | BatFiles, range_: Decimal,
                 cat_only: bool | None = None, workers: int = 1) -> None:
        self.files = files
        self.range_ = range_
        self.cat_only = cat_only
        self.workers = workers

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

        alignment = AlignmentParser(self.files.alignment, self.range_, self.cat_only).parse()
        taxonomy = tax.Taxonomy.from_nodes_file(self.files.nodes)
        fastaid2LCAtaxid = read_fastaid2LCAtaxid(
            self.files.fastaid2LCA, set(alignment.hit_names), workers=self.workers)
        fastaid2taxid = taxonomy.taxa_of_hits(alignment.hit_names, fastaid2LCAtaxid)
        branches = read_branch_taxids(self.files.branches)

        return ClassificationInputs(
            entity_type=entity_type,
            entity2ORFs=entity2ORFs,
            contig_names=contig_names,
            contig2ORFs=contig2ORFs,
            alignment=alignment,
            taxonomy=taxonomy,
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
