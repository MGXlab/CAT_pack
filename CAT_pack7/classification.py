import logging
from time import sleep

from .utils.logging import Status


def contig_classification(args, files, report):
    """placeholder for the classify_contigs function"""
    log = logging.getLogger("CAT_pack")
    for i in range(100):
        log.info(f"Classifying contig {i}")
        report("Classify", Status.RUNNING, i, 100)
        sleep(0.1)
