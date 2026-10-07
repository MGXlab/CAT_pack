import logging
import multiprocessing
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from ..config.settings import BatFiles
from ..io.parsers import FastaParser
from ..utils.errors import InputError
from ..utils.logging import Status
from ..utils.memory import get_memory_budget

log = logging.getLogger("CAT_pack")

def run_protein_prediction(settings, report, predictor):
    if predictor == "pyrodigal":
        return run_pyrodigal(settings, settings.files, report)
    return None


def read_contig_batch(input_paths, max_memory_bytes):
    """Add whole contigs in order.

    Yield before adding a contig that exceeds the estimated memory allowance.
    """
    if max_memory_bytes <= 0:
        raise InputError("The batch memory allowance must be positive.")
    contig_names, sequences, n_bytes = [], [], 0
    for path in input_paths:
        log.info(f"Parsing fasta file: {path}")
        for record in FastaParser(path):
            sequence = record.sequence.encode()
            # The 2048 is an estimation for overhead per contig
            estimated = len(sequence) * 128 + 2048
            if sequences and n_bytes + estimated > max_memory_bytes:
                yield contig_names, sequences
                contig_names, sequences, n_bytes = [], [], 0
            contig_names.append(record.name)
            sequences.append(sequence)
            n_bytes += estimated

    # a catch for the end batch if not filled up completely
    if sequences:
        yield contig_names, sequences


def find_genes(sequence):
    import pyrodigal
    gene_finder = pyrodigal.GeneFinder(meta=True)
    return gene_finder.find_genes(sequence)


def predict_contig(item):
    header, sequence = item
    prediction = find_genes(sequence)
    proteins, gff = StringIO(), StringIO()
    prediction.write_translations(proteins, sequence_id=header)
    prediction.write_gff(gff, sequence_id=header, header=False)
    return proteins.getvalue(), gff.getvalue()


def run_pyrodigal(settings, files, report):
    """Predict proteins from the CAT contigs or the BAT bin files."""

    log.info(f"Running Pyrodigal for ORF prediction. Files {files.proteins_fasta}"
        f" and {files.proteins_gff} will be generated.")

    # TODO: find a nicer solution for this
    input_paths = files.bin_paths if isinstance(files, BatFiles) else (files.contigs,)
    n_contigs = 0
    budget = get_memory_budget(settings.memory, settings.threads)

    # work within the temp directory and after completion copy over to the
    # finished file/folder
    with TemporaryDirectory(prefix=".pyrodigal-prot-prediction-", dir=files.proteins_fasta.parent) as tmpdir:
        proteins_fasta = Path(tmpdir) / files.proteins_fasta.name
        proteins_gff = Path(tmpdir) / files.proteins_gff.name

        with multiprocessing.Pool(processes=budget.workers) as pool:
            with open(proteins_fasta, "w") as outf1, open(proteins_gff, "w") as outf2:
                # setup own headers :)
                # https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
                outf2.write("##gff-version 3\n") # Version must be added following specifications
                outf2.write("# CAT_pack protein prediction; predictor=pyrodigal\n")
                for contig_names, sequences in read_contig_batch(input_paths, budget.batch_bytes):

                    predictions = pool.map(predict_contig, zip(contig_names, sequences), chunksize=1)
                    for protein_text, gff_text in predictions:
                        outf1.write(protein_text)
                        outf2.write(gff_text)

                    n_contigs += len(contig_names)
                    report("Protein prediction", Status.RUNNING, n_contigs, None)
                    # explicid release, to prevent memory growth
                    del predictions, protein_text, gff_text, contig_names, sequences

        proteins_fasta.replace(files.proteins_fasta)
        proteins_gff.replace(files.proteins_gff)
