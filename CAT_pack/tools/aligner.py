"""protein aligners DIAMOND and MMseqs2"""
import logging
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

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


def run_mmseqs2(mmseqs: MMseqsSettings, report):
    log.warning("HIII, as of today I don't know how to run MMseqs2. "
                "Please have a chocolate chip cookie and a good cup of coffee."
                " The porting is in progress.\n")

    log.warning(
    r"""
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


def run_aligner(arguments: DiamondSettings | MMseqsSettings, report, *, _=None):
    if type(arguments) == DiamondSettings:
        log.info(f"Aligning {arguments.query} to {arguments.alignment} with DIAMOND")
        run_diamond(arguments, report)
    elif type(arguments) == MMseqsSettings:
        run_mmseqs2(arguments, report)
        pass
