.. _faq:

##########################
Frequently Asked Questions
##########################

What is an appropriate ``-L/--lnlike-thresh`` value to use?
===================================================
The parameter ``-L`` or ``--lnlike-thresh`` specifies a log-likelihood differential cutoff.  When the three possible topologies for a given quartet are evaluated using molecular data for that quartet by RAxML, three separate likelihood values are generated, which are then log-transformed.  When the cutoff value is specified this means that you want to ensure that the tree with the best likelihood exceeds the next-most likely tree's likelihood by at least the value of the cutoff.  This prevents the method from selecting from two or three trees with nearly indistinguishable likelihoods.  These quartet replicates that pass the likelihood cutoff are tabulated separated and calculated as the QI score.  

On a practical level, increasing the ``-L`` parameter will decrease QI and increase the number of replicates counted as part of QC and QD.  Lowering the ``-L`` parameter will increase the QI score, but increase the number of potentially arbitrary counts in QC and QD.

More statistically, the log-likelihood threshold is equivalent to specifying a minimum likelihood ratio (subtracting log-likelihoods is equivalent to dividing raw likelihoods). Let the null hypothesis be zero difference between the tree likelihoods (i.e., any difference is random noise). For each of the two likelihood estimates compared, if error is distributed in a standard normal distribution N(0,1), then 95% of the distribution lies within ~2 standard deviations. This means the difference in the two log-likelihoods is chi-square distributed with the critical value for 1 degree of freedom at 3.8414 for a 95% confidence cutoff.  Therefore, log-likelihood values that differ by 3.84 divided by 2 (for a two-tailed test) equals 3.84/2 = 1.92, which is approximately equal to 2.  Therefore, a ``-L/--lnlikethresh`` value of 2 conservatively evaluates that the two likelihoods differ significantly at the alpha =0.05 level.
