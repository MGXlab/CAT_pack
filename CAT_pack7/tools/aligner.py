"""aligner script"""
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class MMseqsArgs:
    pass

@dataclass
class DiamondArgs:
    diamond: Path
    mode: str
    query: str
    database: str
    top: int
    no_self_hits: bool
    threads: int
    block_size: int
    index_chunks: int
    tmpdir: str
    compression: bool
    alignment: str
    verbose: bool
    blast_flavour: str = "blastp"
    matrix: str = "BLOSUM62"
    evalue: str = "0.001"

    def get_command(self):
        arguments = [
            self.diamond, self.blast_flavour,
            "-q", self.query,
            "-d", self.database,
            "--top", self.top,
            "--matrix", self.matrix,
            "--evalue", self.evalue,
            "-o", self.alignment,
            "-p", self.threads,
            "--block-size", self.block_size,
            "--index-chunks", self.index_chunks,
            "tmpdir", self.tmpdir,
            "--compress", int(self.compression), # boolean to int forces casting of 1 or 0 and then convert to string, else 'True' will be created
            f"--{self.mode}" if self.mode != 'default' else "",
            "--quiet" if not self.verbose else "",
            "--no-self-hits" if not self.no_self_hits else ""
        ]
        return [str(arg) for arg in arguments] # make sure everything is a string

    # def get_table(self):
    #     for item in self.__annotations__:


def run_diamond(diamond: DiamondArgs, log, report):

    log.write(
        "Homology search with DIAMOND is starting. Please be patient. Do not "
        "forget to cite DIAMOND when using CAT or BAT in your publication.\n"
    )

    # log.write(
    #     f"  blast flavour:     {args.blastp}"
    # )

    # blast_settings = Table(title="BLAST settings")
    # blast_settings.add_column("Setting")
    # blast_settings.add_column("Value")
    # blast_settings.add_row("Blast flavour", args.blastp)
    # blast_settings.add_row("Mode" , args.diamond_mode)
    # blast_settings.add_row("BLAST threads", args.threads)

    if not os.path.isdir(diamond.tmpdir):
        os.mkdir(diamond.tmpdir)
    log.flush()

    try:
        subprocess.check_call(diamond.get_command())
        log.write(f"DIAMOND finished. Alignment written to {diamond.alignment}\n")
        log.flush()
    finally:
        # removal of tmp directory
        shutil.rmtree(diamond.tmpdir, ignore_errors=True)





def run_aligner(arguments: DiamondArgs | MMseqsArgs, log, report):
    if type(arguments) == DiamondArgs:
        run_diamond(arguments, log, report)
    elif type(arguments) == MMseqsArgs:
        #run_mmseqs2(arguments)
        pass