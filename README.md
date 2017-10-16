# CAT

## Introduction

CAT is a pipeline for taxonomic classification of long sequences implemented in Python. It uses the software Prodigal for gene prediction and Diamond as an alignment tool. CAT also requires a database of reference sequences which have NCBI accession numbers in the headers (as contemporary NCBI databases have) and corresponding NCBI taxonomy tree files. CAT can be run from two intermediate steps if files are formated appropriately (see examples).

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

You will need to specify absolute paths to Prodigal and Diamond as well as to the directory with taxonomy tree files inside CAT:

	diamond = '/absolute/path/to/directory/with/diamond/'

	prodigal = '/absolute/path/to/directory/with/prodigal/'

	path_to_taxonomy_files=’/absolute/path/to/files/directory/’

If you made CAT executable and added it in your PATH environment variable, it can be run using the command like this:

	$ CAT --fna sequences.fna --db reference_database.dmnd --prefix library_one

To get help:

	$ CAT -h

## What if I want to use a different ORF caller or aligning tool?

Actually, not a problem. CAT is a universal classification tool and you can apply it to sequences of any origin.

CAT has just a few constrains:

1. Headers of predicted proteins should have the following form:

	 \>NameOfContigWithoutSpaces(symbols like **. : _** are allowed)_**_OrfNumber(only digits, separated from name by one underline)**

2. A file with results of an alignment should be tab-delimited, without an header. Protein name should be in the first column and bitscore value in the last one:

	NameOfContigWithoutSpaces_**_OrfNumber** \<tab\> field2 \<tab\> ... \<tab\> fieldN \<tab\> **Bitscore**


For more details about the algorithm, please see http://biorxiv.org/content/early/2016/09/01/072868

