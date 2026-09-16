"""protein aligners DIAMOND and MMseqs2"""
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

log = logging.getLogger("CAT_pack")

AlignerName = Literal["diamond", "mmseqs2"] #makes 100% sure one of the two is selected


# So the idea is, Arguments are the supplied commando's by the user,
# a dataclass will store them and will always have defaults so there will be
# no unexpected errors there. Then within verification.py the arguments are
# checked and will be set up into a SomeFunctionSettings dataclass. This
# dataclass will then exclusively be used inside CAT_pack.
@dataclass(frozen=True)
class DiamondArgs:
    mode: str = "default"
    no_self_hits: bool = False
    block_size: float = 12.0
    index_chunks: int = 1
    path_to_diamond: Path | None = None


@dataclass(frozen=True)
class MMseqsArgs:
    sensitivity: float = 5.7
    split_memory_limit: str = "0"
    executable: Path | None = None

@dataclass(frozen=True)
class DiamondSettings:
    diamond: Path
    query: Path
    database: Path
    alignment: Path
    tmpdir: Path
    threads: int
    top: int
    mode: str
    no_self_hits: bool
    block_size: float
    index_chunks: int
    compression: bool
    verbose: bool
    blast_flavour: str = "blastp"
    matrix: str = "BLOSUM62"
    evalue: str = "0.001"

    def get_command(self):
        command = [
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
            "--compress", int(self.compression),
            f"--{self.mode}" if self.mode != 'default' else "",
            "--quiet" if not self.verbose else "",
            "--no-self-hits" if not self.no_self_hits else ""
        ]
        return [str(arg) for arg in command] # make sure everything is a string


@dataclass(frozen=True)
class MMseqsSettings:
    pass





def run_diamond(diamond: DiamondSettings, report):
    log.warning(
        "Homology search with DIAMOND is starting. Please be patient. Do not "
        "forget to cite DIAMOND when using CAT or BAT in your publication.\n"
    )

    if not diamond.tmpdir.is_dir():
        log.info(f"making tmp dir: {diamond.tmpdir}")
        diamond.tmpdir.mkdir(parents=True, exist_ok=True)

    log.info(f"Running command: {' '.join(diamond.get_command())}")

    try:
        subprocess.check_call(diamond.get_command())
        log.info(f"DIAMOND finished. Alignment written to {diamond.alignment}\n")
    finally:
        shutil.rmtree(diamond.tmpdir, ignore_errors=True)


def run_mmseqs(mmseqs: MMseqsArgs, report):
    log.warning("HIII, as of today I don't know how to run MMseqs2. "
                "Please have a chocolate chip cookie and a good cup of coffee."
                " The porting is in progress.\n")

    log.warning(
    f"""
               _                        
           \`*-.                    
            )  _`-.                 
           .  : `. .                
           : _   '  \               
           ; *` _.   `*-._          
           `-.-'          `-.       
             ;       `       `.     
             :.       .        \    
             . \  .   :   .-'   .   
             '  `+.;  ;  '      :   
             :  '  |    ;       ;-. 
             ; '   : :`-:     _.`* ;
    \[bug] .*' /  .*' ; .*`- +'  `*' 
          `*-*   `*-*  `*-*'            By: B. Kozlowski
    """)


def run_aligner(arguments: DiamondSettings | MMseqsSettings, report):
    if type(arguments) == DiamondSettings:
        log.info(f"Aligning {arguments.query} to {arguments.alignment} with DIAMOND")
        run_diamond(arguments, report)
    elif type(arguments) == MMseqsSettings:
        #run_mmseqs2(arguments)
        pass