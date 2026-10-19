# def run_bwa_mem(
#         path_to_bwa,
#         path_to_samtools,
#         contigs_fasta,
#         read_file,
#         out_prefix,
#         nproc,
#         log_file
#
import logging
from pathlib import Path

from CAT_pack7.config.settings import RatSettings
from CAT_pack7.tools.tool_runner import run_tool


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
    with tempfile.TemporaryDirectory(prefix="RAT-", dir=files.bam_files[0].parent):
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


        try:
            log.info("Running bwa mem...")

            bwa_command = [
                str(settings.bwa), "mem",
                "-t", str(settings.threads),
                str(files.contigs_fasta),
                str(files.read_files[0])
            ]

            if len(files.read_files) == 2:
                bwa_command.append(str(files.read_files[1]))

            bwa_process = subprocess.Popen(bwa_command, stdout=subprocess.PIPE)


            try:
                # Run Samtools view and remove bwamem file
                command = [
                    str(settings.samtools), "sort",
                    "-@", str(settings.threads/2),
                    "-O", "BAM",
                    "-T", str(settings.tmp / "sort"),
                    "-o", str(partial),
                    "-"
                ]

                samtools_process = subprocess.Popen(command, stdin=bwa_process.stdout)
            finally:
                bwa_process.stdout.close()



                    message = "Samtools view finished abnormally."

                message = "Bwa mem finished abnormally."


            message = "Bwa index finished abnormally."


    message = "Read mapping done!\n"
    return
