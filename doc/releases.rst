.. _releases:

########
Releases
########

v1.2 - Major Update
======
* In order to make all four QS components have a "perfect" score of 1, the Quartet Uncertainty (QU) and Quartet Differential (QD) scores have been inverted. QU is now known Quartet as Informativeness (QI) is the inverse measure (where QI = 1 - QU).  QD is also now inverted in scale, so that 1.0 means no differential in the discordant trees, and 0 means all one discdorant option.
* This version is consistent with the release of the updated bioRxiv manuscript.

v1.1.1
======
* Changed from -A/--amino-acids to -d/--data-type with 'nuc', 'amino', and 'cat' for binary/categorical data
* Amino acid mode should now work in PAUP as well.
* Fixed some issues with temp-dir and results-dir that were causing paths to revert to the root
* **the -e/--temp-dir parameter folder name now MUST contain 'QuartetSampling' to prevent file deletion**

v1.1
====
Contains small fixes to enable efficient documentation and minor changes to argument documentation.

v1.0
====
First public release, concurrent with bioRxiv publication.

