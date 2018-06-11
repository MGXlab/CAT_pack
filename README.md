# CAT

## Introduction

CAT is a pipeline for taxonomic classification of long sequences or bins of sequences implemented in Python. It uses the software Prodigal for gene prediction or translation in 6 frames for sequences with high error rates. Diamond software is used for homology search. CAT also requires a database of reference sequences which have NCBI accession numbers in the headers (as contemporary NCBI databases have) and corresponding NCBI taxonomy tree files. CAT can be run from two intermediate steps if files are formated appropriately (see examples).

## Dependencies and where to get them

Python 3.5.2

diamond	http://github.com/bbuchfink/diamond  (tested on version v0.8.31)

prodigal	http://github.com/hyattpd/Prodigal  (tested on version v2.6.3)

NCBI taxonomy tree files:

From ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ :

from taxdump archive:

 		 names.dmp 
	
 		 nodes.dmp
	
From ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/

 		 prot.accession2taxid

## Getting started

Before you get started, please check that reference sequences have headers in the following format (accession number with version):

>\>WP_003131952.1 <and anything else>

As Diamond saves into output alignment file only the part of the header of a reference sequence that comes before the first space, to preserve information about functional annotation we suggest to change NCBI fasta headers with the following command:

	$ perl -wple 's/^(>\S+\.\d+)\s(.+?) \[.*/$1:$2/g; s/ /_/g; if (length$_ > 200) {s/^(.{200}).*/$1/g};' database > database_formatted

This will change headers from this variant:

>\>WP_003131952.1 30S ribosomal protein S18 [Lactococcus lactis]

to this one:

>\>WP_003131952.1:30S_ribosomal_protein_S18

This is an optional step, CAT will work with the original NCBI headers as well.

After that you can generate a Diamond database from reference sequences as described here http://ab.inf.uni-tuebingen.de/data/software/diamond/download/public/manual.pdf

For convenience, you can specify paths to Prodigal and Diamond as well as to the directory with taxonomy tree files and database inside CAT.

## Usage

If you made CAT executable, added it in your PATH environment variable as well as Diamond and Prodigal, and taxonomy files are in a working directory, it can be run in 'contigs annotation' mode as:

	$ CAT --fna sequences.fna --db reference_database.dmnd --prefix CAT_out

And in 'bins annotation' mode:

	$ CAT --bins /path/to/bins/directory/ --db reference_database.dmnd --prefix CAT_out

To get help:

	$ CAT -h

## What if I want to use a different ORF caller or aligning tool?

Actually, not a problem. CAT is a universal classification tool and you can apply it to sequences of any origin.

If you would like to use protein sequences you already have, you will need to make a tab-delimited file connecting names of contigs (first field) to protein names (second field, one row per protein). 

If you would like to use your own alignment file, please check, that it is a tab-delimited file with information about Bitscore in the last field/column:

	ProteinName \<tab\> field2 \<tab\> ... \<tab\> fieldN \<tab\> **Bitscore**


For more details about the algorithm, please see http://biorxiv.org/content/early/2016/09/01/072868

