"""protein aligners DIAMOND and MMseqs2"""
import logging
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from CAT_pack.io.parsers import optionally_compressed_handle
from CAT_pack.io.writers import top_tabular_alignment
from .tool_runner import run_tool
from ..config.settings import DiamondSettings, MMseqsSettings

log = logging.getLogger("CAT_pack")


def run_diamond(diamond: DiamondSettings, report):
    log.info("Homology search with DIAMOND is starting. Please be patient.")

    tmp_parent = diamond.tmpdir.resolve()
    tmp_parent.mkdir(parents=True, exist_ok=True) # DEVNOTE: Maybe do this already in validation

    with TemporaryDirectory(
        prefix="CAT_pack_DIAMOND_", dir=tmp_parent, ignore_cleanup_errors=True,
    ) as tmpdir:
        diamond = replace(diamond, tmpdir=Path(tmpdir))
        run_tool(diamond.get_command(), "DIAMOND", operation="homology search")
        log.info(f"DIAMOND finished. Alignment written to {diamond.alignment}\n")

def run_mmseqs2(mmseqs: MMseqsSettings, report) -> None:
    log.info("Homology search with MMseqs2 is starting. Please be patient.")
    if mmseqs.gpu and mmseqs.sensitivity:
        log.warning("MMseqs2 GPU search uses maximum sensitivity; --sensitivity does not apply.")

    tmp_parent = mmseqs.tmpdir.resolve()
    tmp_parent.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(
        prefix="CAT_pack_MMseqs2_", dir=tmp_parent, ignore_cleanup_errors=True,
    ) as tmpdir:
        tmpdir = Path(tmpdir)
        tmp_alignment_file = tmpdir / "alignment.tsv"
        run_settings = replace(mmseqs, alignment=tmp_alignment_file, tmpdir=tmpdir)
        run_tool(run_settings.get_command(), "MMseqs2", operation="homology search")

        log.info("Subsetting alignment table based on --top %s.", mmseqs.top)

        with (
            optionally_compressed_handle(tmp_alignment_file, 'r') as f1,
            optionally_compressed_handle(mmseqs.alignment, "wt") as outf1,
        ):
            top_tabular_alignment(
                f1, outf1, mmseqs.top, 2 if mmseqs.tiny_out else 11,
            )

    log.info("MMseqs2 finished. Alignment written to %s\n", mmseqs.alignment)


def run_aligner(arguments: DiamondSettings | MMseqsSettings, report):
    if type(arguments) == DiamondSettings:
        log.info(f"Aligning {arguments.query} to {arguments.alignment} with DIAMOND")
        run_diamond(arguments, report)
    elif type(arguments) == MMseqsSettings:
        run_mmseqs2(arguments, report)
