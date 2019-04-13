# Changelog

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
