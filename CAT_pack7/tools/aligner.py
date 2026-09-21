"""protein aligners DIAMOND and MMseqs2"""
import logging
import shutil
import subprocess

from ..settings import DiamondSettings, MMseqsSettings

log = logging.getLogger("CAT_pack")


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


def run_mmseqs2(mmseqs: MMseqsSettings, report):
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
        run_mmseqs2(arguments, report)
        pass
