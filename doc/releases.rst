.. _releases:

########
Releases
########

v1.3.1 - Added RESULTS.node.counts.csv
======================================
The output now includes a requested feature where the counts of all three quartet options are shown in a separate file alongside representative quartet topologies with labels.  In cases where QD is high, this will allow you to check which discordant option is more common by looking at the counts.

v1.3 - Change to RAxML-ng, modernization of flags
=================================================
* **The default is now RAxML-ng**, though classic RAxML 8.1+ is still available.  Performance of QS in RAxML-ng on test datasets gives similar results, but is much more efficient.
* **Discontinuation of single-letter flags.** Single-letter flags are no longer functional, but have been retained in the comments of the main *quartet_sampling.py* script to allow forward updating.  This change brings QS into line with evolving commnuity best practices that find better clarity and reproducibility when word-based flags are used exclusively.
* **--engine syntax change:** As we expand support for new tree inference and likelihood evaluation engines, this flag is a more flexible alternative to the previous named ones.  See also notes on *--engine-exec* and *--engine-model*.
* Support for IQ-TREE is also available FOR TESTING PURPOSES ONLY!  We have not fully evaluated QS under the IQ-TREE engine, so use this with caution and be sure to document this in your methods.

v1.2.1 - Efficiency update
==========================
* RAxML now evaluates all three tree configurations in the same program call.  This reduces RAxML calls and increases efficiency.
* Clarification in the manual that RAxML 8.1+ is required for the program to execute correctly.

v1.2 - Major Update
===================
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

