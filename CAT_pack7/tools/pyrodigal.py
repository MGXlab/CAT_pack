import logging
import multiprocessing

from ..parsers import FastaParser
from ..settings import BatFiles


def run_protein_prediction(settings, report, predictor):
    if predictor == "pyrodigal":
        return run_pyrodigal(settings, settings.files, report)
    return None


def run_pyrodigal(settings, files, report):
    """Predict proteins from the CAT contigs or the BAT bin files."""
    log = logging.getLogger("CAT_pack")
    import pyrodigal
    gene_finder = pyrodigal.GeneFinder(meta=True)

    log.warning(f"Running Pyrodigal for ORF prediction. Files {files.proteins_fasta}"
        f" and {files.proteins_gff} will be generated. Do not forget to cite"
        " Pyrodigal and Prodigal when using CAT or BAT in your publication.")

    input_paths = files.bin_paths if isinstance(files, BatFiles) else (files.contigs,)
    contig_names = []
    sequences = []
    for path in input_paths:
        log.info(f"Parsing contigs fasta {path}")
        for record in FastaParser(path):
            contig_names.append(record.name)
            sequences.append(record.sequence.encode())

    with multiprocessing.Pool(processes=settings.threads) as pool:
        predictions = pool.map(
            gene_finder.find_genes,
            sequences,
        )

    with open(files.proteins_fasta, "w") as outf1, open(files.proteins_gff, "w") as outf2:
        for header, prediction in zip(contig_names, predictions):
            prediction.write_translations(outf1, sequence_id=header)
            prediction.write_gff(outf2, sequence_id=header)
