import logging
import multiprocessing
from pathlib import Path
from tempfile import TemporaryDirectory

from ..parsers import FastaParser
from ..settings import BatFiles
from ..utils.logging import Status

log = logging.getLogger("CAT_pack")

def run_protein_prediction(settings, report, predictor):
    if predictor == "pyrodigal":
        return run_pyrodigal(settings, settings.files, report)
    return None


def read_contig_batch(input_paths, max_contigs, max_bases):
    """Add whole contigs in order, yield when the batch amound or size limit is reached

    It's allowed to go a bit over the size limit, this can be made into a hard limit
    but keeping in mind that diamond does use a lot of memory, this doesn't impact it that much.
    """
    contig_names, sequences, n_bases = [], [], 0
    for path in input_paths:
        log.info(f"Parsing fasta file: {path}")
        for record in FastaParser(path):
            sequence = record.sequence.encode()
            contig_names.append(record.name)
            sequences.append(sequence)
            n_bases += len(sequence)
            if len(sequences) >= max_contigs or n_bases >= max_bases:
                yield contig_names, sequences
                contig_names, sequences, n_bases = [], [], 0 # Reset :)

    # a catch for the end batch if not filled up completely
    if sequences:
        yield contig_names, sequences


def find_genes(sequence):
    import pyrodigal
    gene_finder = pyrodigal.GeneFinder(meta=True)
    return gene_finder.find_genes(sequence)


def run_pyrodigal(settings, files, report):
    """Predict proteins from the CAT contigs or the BAT bin files."""

    log.info(f"Running Pyrodigal for ORF prediction. Files {files.proteins_fasta}"
        f" and {files.proteins_gff} will be generated.")

    # TODO: find a nicer solution for this
    input_paths = files.bin_paths if isinstance(files, BatFiles) else (files.contigs,)
    n_contigs = 0

    # work within the temp directory and after completion copy over to the
    # finished file/folder
    with TemporaryDirectory(prefix=".pyrodigal-prot-prediction-", dir=files.proteins_fasta.parent) as tmpdir:
        proteins_fasta = Path(tmpdir) / files.proteins_fasta.name
        proteins_gff = Path(tmpdir) / files.proteins_gff.name

        with multiprocessing.Pool(processes=settings.threads) as pool:
            with open(proteins_fasta, "w") as outf1, open(proteins_gff, "w") as outf2:
                # setup own headers :)
                # https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
                outf2.write("##gff-version 3\n") # Version must be added following specifications
                outf2.write("# CAT_pack protein prediction; predictor=pyrodigal\n")
                for contig_names, sequences in read_contig_batch(
                    input_paths, settings.pyrodigal.max_n_contigs_p_batch,
                    settings.pyrodigal.max_n_bases_p_batch):

                    predictions = pool.map(find_genes, sequences)
                    for header, prediction in zip(contig_names, predictions):
                        prediction.write_translations(outf1, sequence_id=header)
                        prediction.write_gff(outf2, sequence_id=header, header=False)

                    n_contigs += len(contig_names)
                    report("Protein prediction", Status.RUNNING, n_contigs, None)
                    # explicid release, to prevent memory growth
                    del predictions, prediction, contig_names, sequences

        proteins_fasta.replace(files.proteins_fasta)
        proteins_gff.replace(files.proteins_gff)
