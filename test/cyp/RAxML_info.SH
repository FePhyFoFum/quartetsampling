

WARNING: RAxML is not checking sequences for duplicate seqs and sites with missing data!



Using BFGS method to optimize GTR rate parameters, to disable this specify "--no-bfgs" 


Alignment has 2379 completely undetermined sites that will be automatically removed from the input data



This is RAxML version 8.1.21 released by Alexandros Stamatakis on May 21 2015.

With greatly appreciated code contributions by:
Andre Aberer      (HITS)
Simon Berger      (HITS)
Alexey Kozlov     (HITS)
Kassian Kobert    (HITS)
David Dao         (KIT and HITS)
Nick Pattengale   (Sandia)
Wayne Pfeiffer    (SDSC)
Akifumi S. Tanabe (NRIFS)

Alignment has 7712 distinct alignment patterns

Proportion of gaps and completely undetermined characters in this alignment: 81.49%

RAxML computation of SH-like support values on a given tree

Using 1 distinct models/data partitions with joint branch length optimization


All free model parameters will be estimated by RAxML
GAMMA model of rate heteorgeneity, ML estimate of alpha-parameter

GAMMA Model parameters will be estimated up to an accuracy of 0.1000000000 Log Likelihood units

Partition: 0
Alignment Patterns: 7712
Name: No Name Provided
DataType: DNA
Substitution Matrix: GTR




RAxML was called as follows:

raxml -f J -m GTRGAMMA -t cyp.tre -n SH -p 12345 -s cyp.phy 


Time after model optimization: 107.164181
Initial Likelihood -168987.306517

NNI interchanges 13 Likelihood -168980.633091
NNI interchanges 1 Likelihood -168980.345819
NNI interchanges 0 Likelihood -168980.337178

Final Likelihood of NNI-optimized tree: -168980.337178

RAxML NNI-optimized tree written to file: /home/smitty/Dropbox/projects/taxonjackknife/quartet_sampling/test/test_cyp/RAxML_fastTree.SH

Same tree with SH-like supports written to file: /home/smitty/Dropbox/projects/taxonjackknife/quartet_sampling/test/test_cyp/RAxML_fastTreeSH_Support.SH

Total execution time: 169.953851
