""""IO (input/output) module of CAT pack

This module consist of four files `download.py`. `parsers.py`, `reference_mapping.py`, `writers.py`

download.py
NOTE: This python script isn't optimized yet and does not follow the parsers.py, writers.py logic yet. 
This python script handles the downloading and processing of the database files

parsers.py
This file handles the parsing of numerous different file types and folders: 
 - FASTA and FASTQ: Can be either nonziped or .gz or .bz2
 - BinParser: Parsing of bins
 - AlignmentParser: DIOMOND aligment parser 
 - ClassificationParser: Load ORF groups, accepted alignment hits, and their taxonomy data
"""