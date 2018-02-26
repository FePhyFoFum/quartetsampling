.. _intro:

###############
Getting Started
###############

What is Quartet Sampling?
=========================
Quartet Sampling (QS) is a method for quantifying branch support values for a phylogenetic tree.  The software requires an input tree topology (branch lengths not required or used) and a molecular alignment in phylip format.  A RAxML-formatted partition file may (optionally) be also used.  The Quartet Sampling method takes a phylogenetic topology and considers each internal (i.e., non-terminal) branch.  Each internal branch has four branches connected to it that lead to one or more taxa.  For each replicate, Quartet Sampling randomly selects one taxon from each of the four groups and evaluates the likelihoods of the three possible topological arrangements of the four taxa.  These replicates then form the basis of branch scores for the QS method.  Full details about the scientific background and equations behind QS can be found in the Quartet Sampling paper at <http://biorxiv.org/content/early/2017/06/10/148536.1>.

How do I cite this program?
===========================
James B Pease, Joseph W Brown, Joseph F Walker, Cody E Hinchliff, Stephen A Smith. 2018. Quartet Sampling distinguishes lack of support from conflicting support in the green plant tree of life. American Journal of Botany. doi:10.1002/ajb2.1016

Please also include the URL <https://www.github.com/fephyfofum/quartetsampling> in your methods section where the program is referenced.

Installation
============
No installation is required, quartetsampling scripts should work as long as Python and RAxML are installed.  The repository can be cloned or downloaded as a .zip file from GitHub.

Requirements
------------
* Python 3.x (2.7 should also work, but 3.x recommended) https://www.python.org/downloads/
* RAxML 8.x (7.x should also work, but 8.x recommended) https://sco.h-its.org/exelixis/web/software/raxml/index.html

**Optional, but recommended:**

* PAUP  http://people.sc.fsu.edu/~dswofford/paup_test/ (required only if PAUP is used instead of RAxML)
* Numpy http://www.numpy.org (required only for ``calc_qs_stats.py``)
* Scipy https://www.scipy.org (required only if ``--calc-qdstats.py`` is activated)
* Figtree http://tree.bio.ed.ac.uk/software/figtree/ (required to view FigTree output)

Preparing your data
===================

Phylogeny
---------

A phylogenetic tree in Newick (parenthetical) format should be used.  Branch lengths are optional and will be ignored by the program.  Removal of support scores in square brackets is recommended.  

.. note:: Internal branch labels will be replaced in output trees.

Sequence Alignment
------------------

An alignment in Relaxed Phylip format (such as used for RAxML) is required.  The alignment can be DNA nucleotides or amino acids (use ``-A/--amino-acid`` in that case). If you have an alignment in FASTA format, we also include a script called ``utils/fasta2phy.py`` that will convert FASTA alignments to Relaxed Phylip. 

.. important:: Labels for the alignment sequences must match the labels on tree terminal branches exactly. All tips in the phylogeny must have a sequence represented in the alignment.  Sequences appearing in the alignment, but not in the tree are allowed, and will be ignored.

Partitions File (Optional)
--------------------------

A partitions file in the style of RAxML may also be used either for a partitioned analysis of likelihood or separate gene tree evaluation.  DNA partitions should use the ``DNA`` prefix, and proteins should use the ``WAG`` prefix instead of ``DNA`` in the example below.  Note that partition names are arbitrary in this case and that number ranges are inclusive.

Example file::
  DNA, p1=1-30
  DNA, p2=31-60
  DNA, p3=61-90
  ...

.. important:: If your alignment is sparse and you use partitioned most, such that you will frequently have partitions with no sequence data for randomly selected quartets of taxa, you may wish to invoke the ``--ignore-error`` option in ``quartet_sampling.py`` to ignore the RAxML errors that will result from these empty partitions.

Basic usage
===========

::
  python quartet_sampling.py -t TREE.nwk -a ALIGNMENT.phy -N 100 -T 4 -L 2

Required Paramters
------------------
``-t/--tree``: File containing a single Newick-formatted phylogeny.  

``-a/--alignment`` Phylip-formatted alignment containing one sequence for each tip in the phylogeny provided.

``-N/--number-of-reps`` Number of replicates per branch to perform.

``-T/--number-of-threads`` Number of threads for Python to use in parallel (does not pass through to RAxML, this is for Python multiprocessing of per-branch replicates in parallel)

Recommended Parameters
----------------------

``-L/--lnlike-thresh``: log-likelihood threshold cutoff, determines the minimum difference by which the best likelihood tree must exceed the second-best likelihood tree when comparing the three possible topologies for a given quartet replicate.

.. note:: If ``-L/--lnlike-thresh`` is omitted, this will invoke an alternative mode where a tree is simply inferred from the sequence data by RAxML (or PAUP*) and likelihoods are not evaluated.  This will result in Quartet Informativeness (QI) scores of 'NA' for all branches.




