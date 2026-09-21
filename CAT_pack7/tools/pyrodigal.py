import logging
import multiprocessing


def run_protein_prediction(settings, report, predictor):
    if predictor == "pyrodigal":
        return run_pyrodigal(settings, settings.files, report)
    return None


def run_pyrodigal(settings, files, report):
    """placeholder for the pyrodigal"""
    log = logging.getLogger("CAT_pack")
    import pyrodigal
    gene_finder = pyrodigal.GeneFinder(meta=True)

    log.warning(f"Running Pyrodigal for ORF prediction. Files {files.proteins_fasta}"
        f" and {files.proteins_gff} will be generated. Do not forget to cite"
        " Pyrodigal and Prodigal when using CAT or BAT in your publication.")

    header2seq = {}
    contig_order = []

    def import_contigs(contigs_fasta):
        log.info(f"Parsing contigs fasta {contigs_fasta}")
        with open(contigs_fasta, "r") as f1:
            for line in f1:
                line = line.rstrip()

                if line.startswith(">"):
                    header = line.rstrip().split(" ")[0].lstrip(">")
                    header2seq.setdefault(header, "")
                    contig_order.append(header)
                else:
                    header2seq[header] += line.rstrip()

    import_contigs(files.contigs)

    with multiprocessing.Pool(processes=settings.threads) as pool:
        predictions = pool.map(
            gene_finder.find_genes,
            [bytes(header2seq[header].encode()) for header in contig_order]
        )

    with open(files.proteins_fasta, "w") as outf1, open(files.proteins_gff, "w") as outf2:
        for header, prediction in zip(contig_order, predictions):
            prediction.write_translations(outf1, sequence_id=header)
            prediction.write_gff(outf2, sequence_id=header)
