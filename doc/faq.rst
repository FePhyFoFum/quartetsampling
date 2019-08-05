.. _faq:

##########################
Frequently Asked Questions
##########################

What is an appropriate ``--lnlike`` value to use?
===========================================================
The parameter ``--lnlike`` or ``--lnlike-thresh`` specifies a log-likelihood differential cutoff.  When the three possible topologies for a given quartet are evaluated using molecular data for that quartet by RAxML, three separate likelihood values are generated, which are then log-transformed.  When the cutoff value is specified this means that you want to ensure that the tree with the best likelihood exceeds the next-most likely tree's likelihood by at least the value of the cutoff.  This prevents the method from selecting from two or three trees with nearly indistinguishable likelihoods.  These quartet replicates that pass the likelihood cutoff are tabulated separated and calculated as the QI score.  

On a practical level, increasing the ``--lnlike`` parameter will decrease QI and increase the number of replicates counted as part of QC and QD.  Lowering the ``--lnlike`` parameter will increase the QI score, but increase the number of potentially arbitrary counts in QC and QD.

More statistically, the log-likelihood threshold is equivalent to specifying a minimum likelihood ratio (subtracting log-likelihoods is equivalent to dividing raw likelihoods). Let the null hypothesis be zero difference between the tree likelihoods (i.e., any difference is random noise). For each of the two likelihood estimates compared, if error is distributed in a standard normal distribution N(0,1), then 95% of the distribution lies within ~2 standard deviations. This means the difference in the two log-likelihoods is chi-square distributed with the critical value for 1 degree of freedom at 3.8414 for a 95% confidence cutoff.  Therefore, log-likelihood values that differ by 3.84 divided by 2 (for a two-tailed test) equals 3.84/2 = 1.92, which is approximately equal to 2.  Therefore, a ``--lnlike`` value of 2 conservatively evaluates that the two likelihoods differ significantly at the alpha =0.05 level.


My partitioned data is not working properly, and is a sparse partitioned set.
=============================================================================
RAxML and RAxML-ng both require that none of data partitions being evaluated are empty.  This creates a unique challenge for quartet samplingon large-sparse matrices with partitioned data, since with sparse data and partitions there might be many partitions of a given quartet that contain insufficient data.  We are working on a solution, but for now recommend using gene-tree mode or unpartitioned on sparse data.
