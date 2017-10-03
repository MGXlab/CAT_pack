# CAT

## Introduction

CAT is a pipeline for taxonomic classification of long sequences implemented in Python. It uses Prodigal software for gene prediction and Diamond alignment tool. CAT also requires database of reference sequences which have NCBI accession numbers in headers (contemporary NCBI databases) and NCBI taxonomy tree files. CAT allows to input files at two intermediate steps if files formatted accordingly (see examples of files).

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

Before start, please check that reference sequences have headers in the following format (accession number with version):

>\>WP_003131952.1 <and anything else>

As Diamond save into alignment file only part of header of reference sequence before the first space, to preserve information about functional annotation we can offer to change NCBI fasta headers with the following command:

	$ perl -wple 's/^(>\S+\.\d+)\s(.+?) \[.*/$1:$2/g; s/ /_/g; if (length$_ > 200) {s/^(.{200}).*/$1/g};' database > database_formatted

This will change headers from this variant:

>\>WP_003131952.1 30S ribosomal protein S18 [Lactococcus lactis]

to this one:

>\>WP_003131952.1:30S_ribosomal_protein_S18

After that you can generate Diamond database as described here http://ab.inf.uni-tuebingen.de/data/software/diamond/download/public/manual.pdf

When dependencies will be downloaded and installed, you need to specify absolute paths to Prodigal and Diamond as well as to taxonomy tree files inside CAT:

	diamond = '/absolute/path/to/executables/of/diamond/'

	prodigal = '/absolute/path/to/executables/of/prodigal/'

	path_to_taxonomy_files=’/absolute/path/to/files/directory/’

If you made CAT executable and added into PATH environment variable, it could be run using the command like this:

	$ CAT --fna sequences.fna --db reference_database.dmnd --prefix library_one

To get help:

	$ CAT -h


For more details about analysis algorithm, please see http://biorxiv.org/content/early/2016/09/01/072868

