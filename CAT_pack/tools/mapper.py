import logging
import subprocess
import tempfile
from pathlib import Path

from CAT_pack.config.settings import RatSettings
from CAT_pack.tools.tool_runner import run_tool, tool_output
from CAT_pack.utils.errors import ExternalToolError


def check_index(contigs_fasta: Path, mapper: str = "bwa") -> bool:
    suffixes = ((".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
                if mapper == "bwa-mem2" else (".amb", ".ann", ".bwt", ".pac", ".sa"))
    return all(
        Path(f"{contigs_fasta}{suffix}").is_file()
        for suffix in suffixes
    )




def run_bwa_mem(settings: RatSettings, report)-> None:
    files = settings.files
    tool = "BWA-MEM2" if settings.mapper == "bwa-mem2" else "BWA MEM"

    log = logging.getLogger("CAT_pack")
    log.info(
        "Running %s for read mapping. Output: %s. Please cite %s and SAMtools "
        "when using RAT in your publication!", tool, files.bam_files[0], tool)
    with (
        tempfile.TemporaryDirectory(prefix="RAT-", dir=files.bam_files[0].parent) as tmp,
        tool_output(tool) as bwa_output,
        tool_output("SAMtools") as samtools_output,
    ):
        try:
            if check_index(files.contigs_fasta, settings.mapper):
                log.info("Contigs FASTA is already indexed.")
            else:
                log.info("Indexing contigs FASTA...")
                run_tool([str(settings.bwa), "index", str(files.contigs_fasta)],
                    tool,"index")
        except Exception as error:
            log.exception(error)
            raise ExternalToolError(tool, f"Indexing failed: {error}") from error

        tmp = Path(tmp)
        tmp_bam_file = tmp / "sorted.bam"
        bwa_process = None
        samtools_process = None
        try:
            log.info("Running %s...", tool)

            bwa_command = [
                str(settings.bwa), "mem",
                "-t", str(max(1, settings.threads // 2)),
                str(files.contigs_fasta),
                str(files.read_files[0])
            ]

            if len(files.read_files) == 2:
                bwa_command.append(str(files.read_files[1]))

            bwa_process = subprocess.Popen(bwa_command, stdout=subprocess.PIPE,
                                           stderr=bwa_output)


            try:
                command = [
                    str(settings.samtools), "sort",
                    "-@", str(max(1, settings.threads // 2)),
                    "-O", "BAM",
                    "-T", str(tmp / "sort"),
                    "-o", str(tmp_bam_file),
                    "-"
                ]
                if settings.sort_memory_mib is not None:
                    command[2:2] = ["-m", f"{settings.sort_memory_mib}M"]

                samtools_process = subprocess.Popen(command, stdin=bwa_process.stdout,
                                                    stdout=samtools_output,
                                                    stderr=samtools_output)
            finally:
                bwa_process.stdout.close()

            samtool_status = samtools_process.wait()
            bwa_status = bwa_process.wait()

            if bwa_status != 0:
                raise ExternalToolError(tool, f"Mapping failed with exit code {bwa_status}.")
            if samtool_status != 0:
                raise ExternalToolError("Samtools", f"Mapping failed due to an unexpected error during samtools running. Error {samtool_status}")

            tmp_bam_file.replace(files.bam_files[0])

        except (OSError, subprocess.CalledProcessError) as error:
            raise ExternalToolError(f"{tool}/SAMtools",
                f"Mapping failed with error: {error}")

        finally:
            for process in (samtools_process, bwa_process):
                if process is not None and process.poll() is None:
                    process.kill()

            for process in (samtools_process, bwa_process):
                if process is not None:
                    process.wait()

    return None
