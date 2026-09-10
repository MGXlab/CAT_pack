from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from shared import run_CAT


@dataclass(frozen=True)
class CatArgs:
    #all the required
    contigs: Path
    database: Path
    taxonomy: Path
    # some of the optional args
    _range: Decimal
    output_prefix: Path
    threads: int = 1

def run_cat(args):
    run_CAT(args)
    return
# shared.run_CAT(args, args.contigs_fasta, args.database_folder,
#                            args.taxonomy_folder, args.log_file, args.quiet,
#                            args.nproc, args.f, args.r, args.out_prefix,
#                            path_to_CAT)
#     message = "Running CAT."
#     give_user_feedback(message, log_file, quiet, show_time=True)