import logging
import subprocess
import tempfile
from pathlib import Path

from CAT_pack7.config.settings import RatSettings
from CAT_pack7.tools.tool_runner import run_tool
from CAT_pack7.utils.errors import ExternalToolError


def check_index(contigs_fasta: Path) -> bool:
    return all(
        Path(f"{contigs_fasta}{suffix}").is_file()
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")
    )




def run_bwa_mem(settings: RatSettings, report)-> None:
    files = settings.files

    log = logging.getLogger("CAT_pack")
    log.info(
        "Running bwa mem for read mapping. File {0}.sorted will be "
        "generated. Do not forget to cite bwa mem and Samtools when using "
        "RAT in your publication!".format(files.mapping_bam))
    with tempfile.TemporaryDirectory(prefix="RAT-", dir=files.bam_files[0].parent) as tmp:
        try:
            if check_index(files.contigs_fasta):
                log.info("Contigs FASTA is already indexed.")
            else:
                log.info("Indexing contigs FASTA...")
                run_tool([str(settings.bwa), "index", str(files.contigs_fasta)],
                    "BWA","index")
        except Exception as error:
            log.exception(error)
            raise ExternalToolError("BWA MEM", "An error occured during indexing (or checking for indexing): {error}")

        tmp = Path(tmp)
        tmp_bam_file = tmp / "sorted.bam"
        try:
            log.info("Running bwa mem...")

            bwa_command = [
                str(settings.bwa), "mem",
                "-t", str(settings.threads/2),
                str(files.contigs_fasta),
                str(files.read_files[0])
            ]

            if len(files.read_files) == 2:
                bwa_command.append(str(files.read_files[1]))

            bwa_process = subprocess.Popen(bwa_command, stdout=subprocess.PIPE)


            try:
                command = [
                    str(settings.samtools), "sort",
                    "-@", str(settings.threads/2),
                    "-O", "BAM",
                    "-T", str(settings.tmp / "sort"),
                    "-o", str(tmp_bam_file),
                    "-"
                ]

                samtools_process = subprocess.Popen(command, stdin=bwa_process.stdout)
            finally:
                bwa_process.stdout.close()

            samtool_status = samtools_process.wait()
            bwa_status = bwa_process.wait()

            if bwa_status != 0:
                raise ExternalToolError("BWA MEM", f"Mapping failed due to an unexpected error during bwa mem running. Error: {bwa_status}")
            if samtool_status != 0:
                raise ExternalToolError("Samtools", f"Mapping failed due to an unexpected error during samtools running. Error {samtool_status}")

            tmp_bam_file.replace(files.bam_files[0])

        except (OSError, subprocess.CalledProcessError) as error:
            raise ExternalToolError("BWA MEM/SAMtools",
                f"Mapping failed with error: {error}")

        finally:
            for process in (samtools_process, bwa_process):
                if process is not None and process.poll() is None:
                    process.kill()

            for process in (samtools_process, bwa_process):
                if process is not None:
                    process.wait()

    return None