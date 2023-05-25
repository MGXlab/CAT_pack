# Changelog

## 5.X
* GTDB support.
* Sequence databases (NCBI nr or GTDB) can be downloaded with `CAT download`, and CAT databases constructed with `CAT prepare`.
* Sensible defaults of DIAMOND parameters for high memory machines: `--top 11 --block_size 12 --index_chunks 1`.
* Preparations for Read Annotation Tool (RAT).

## 5.2.3
Minor bug fix for `CAT add_names`.

## 5.2.2
We have added the DIAMOND specific `--no_self_hits` flag. We have also added some extra checks and removed redundancy from the parser code. Databases constructed by `CAT prepare` now have a slightly different naming scheme.

## 5.2.1
Minor bug fix for `CAT prepare`.

## 5.2
`CAT prepare` now uses the latest taxonomy mapping files from NCBI, significantly expanding taxonomic coverage of proteins in nr. File integrity of downloads is assessed based on md5 checksums. The ORF2LCA output file contains a new column for the number of hits the classification is based on. We have made textual changes to the output files to better reflect the meaning of 'classified' and 'not classified' in different contexts.

## 5.1.2
Code streamlining.

## 5.1.1
CAT and BAT can now compress the DIAMOND alignment file, and import gzip compressed alignment files.

## 5.1
The code has been rewritten to prepare for future extensions. We have also added the `--verbose` flag.

## 5.0.5.
Skip hidden files in bin folder.

## 5.0.4.
We have added the `--no_stars` flag alongside a minor bug fix.

## 5.0.3
Bug fix for single bin mode.

## 5.0.2
Floating point numbers have been changed to decimals.

## 5.0.1
Updated license to MIT.

## 5.0
We have simplified the output table format: we have added a 'reason' column, which shows the number of ORFs a classification is based on and the total number of predicted ORFs on a contig/MAG. In case of an unclassified sequence, the reason for this is shown in this column as well. Moreover, `add_names` now has an option to exclude the bit-score support scores from the lineage!

## 4.6
We have added the DIAMOND `--top` parameter and the `--I_know_what_Im_doing` flag for experimental features.

## 4.5
BAT can now be run in single bin mode. The familiar `./CAT bins` is still the go-to option if you want to classify multiple MAGs, but if it's only one MAG you are interested in try out `./CAT bin`! An added benefit of single bin mode is that you can use the alignment and predicted protein files of the BAT run to classify individual contigs within the MAG with CAT, or the other way around.

## 4.4
We have added DIAMOND specific options. This allows you to use sensitive mode, and tune memory and temporary disk space usage during alignment! Moreover, you can now force CAT and BAT to overwrite existing files.

## 4.3.4
We extended some of the pre-flight checks.

## 4.3.3
Minor bug fix.

## 4.3.2
A fruity update: CAT and BAT are now macOS compatible!

## 4.3.1
We removed the psutil dependency.

## 4.3
Prepare now checks whether the RAM of your computer is large enough. If not, not to worry! We have put preconstructed databases online.

## 4.2
Code streamlining.

## 4.1
CAT and BAT leave much less footprints: the size of the alignment output files is greatly reduced, alignment is now up to 3 times faster than previous releases.

## 4.0
CAT and BAT have been completely rewritten, bumping the version up to 4.0.
