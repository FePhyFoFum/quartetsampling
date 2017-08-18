Program Parameter Descriptions
##############################

.. quartet_sampling:

quartet_sampling
================

Description
-----------

quartet_samping.py: Quartet Sampling method for
phylogenetic branch support evaluation

<http://www.github.com/FePhyFoFum/quartetsampling>


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-a/--alignment`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Alignment file in "relaxed phylip" format, as used by RAxML.

**Type:** file path; **Default:** None



``-t/--tree`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The input tree in Newick (parenthetical) format.

**Type:** file path; **Default:** None



``-N/--number-of-reps`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The number of replicate quartet topology searches to be performed at each node.

**Type:** integer; **Default:** 100



``-T/--number-of-threads`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The number of parallel threads to be used by Python for quartet topology searches.

**Type:** integer; **Default:** 1



``--calc-qdstats``
^^^^^^^^^^^^^^^^^^

**Description:** EXPERIMENTAL: Calculates Chi-square test for QD tree frequencies. Use only  if Scipy is available. Will increase running time.

**Type:** boolean flag



``-d/--data-type``
^^^^^^^^^^^^^^^^^^

**Description:** (nuc)leotide, (amino) acid, or (cat)egorical data

**Type:** None; **Default:** ['nuc']

**Choices:** ('nuc', 'amino', 'cat')


``-e/--temp-dir``
^^^^^^^^^^^^^^^^^

**Description:** A directory to which temporary files will be saved. If not supplied, 'QuartetSampling' will be created in the current working directory. When specifying a custom temporary output the characters 'QuartetSampling' must appear in the directory name to prevent accidental file deletion. (default='./QuartetSampling'

**Type:** file path; **Default:** None



``-g/--genetrees``
^^^^^^^^^^^^^^^^^^

**Description:** Use partitions file (RAxML format) to divide the alignment into separate gene tree regions. Gene alignments will be sampled random for the quartet topology searches.

**Type:** file path; **Default:** None



``--ignore-errors``
^^^^^^^^^^^^^^^^^^^

**Description:** Ignore RAxML and PAUP erroneous runs

**Type:** boolean flag



``--low-mem``
^^^^^^^^^^^^^

**Description:** Do not store large alignment in memory for whole-alignment (non-genetree) mode

**Type:** boolean flag



``--max-random-sample-proportion``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The proportion of possible replicates explored unsuccessfully by the random generation procedure before it gives up. Because this generates random replicates, it takes progressively longer as it proceeds. To avoid long runtimes, the recommended range is < 0.5 (which is the default).

**Type:** float; **Default:** None



``-o/--results-dir``
^^^^^^^^^^^^^^^^^^^^

**Description:** A directory to which output files will be saved. If not supplied, the current working directory will be used. (default is current folder).

**Type:** file path; **Default:** None



``-p/--stop-node-number``
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** An integer denoting the node at which to stop. Will include nodes with indices <= the stop node number. This argument may be used to limit the length of a given run in case only a certain part of the tree is of interest. Nodes will be read from topologically identical (and isomorphic!) input trees in deterministic order.

**Type:** integer; **Default:** None



``--paup-executable``
^^^^^^^^^^^^^^^^^^^^^

**Description:** The name or path of the PAUP executable to be used for calculated quartets.

**Type:** None; **Default:** ['paup']



``-q/--partitions``
^^^^^^^^^^^^^^^^^^^

**Description:** Partitions file in RAxML format. If omitted then the entire alignment will be treated as one partition for all quartet replicate topology searches.

**Type:** file path; **Default:** None



``-r/--result-prefix``
^^^^^^^^^^^^^^^^^^^^^^

**Description:** A prefix to put on the result files.

**Type:** string; **Default:** None



``--raxml-model``
^^^^^^^^^^^^^^^^^

**Description:** Advanced: specify a custom RAxML model name for the raxml '-m' parameter

**Type:** None; **Default:** None



``--retain-temp``
^^^^^^^^^^^^^^^^^

**Description:** Do not remove temporary files

**Type:** boolean flag



``-s/--start-node-number``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** An integer denoting the node to which to start from. Nodes will be read from topologically identical (and isomorphic!) input trees in deterministic order, so this argument may be  used to restart at an intermediate position (in case the previous run was canceled before completion, for example).

**Type:** integer; **Default:** None



``-v/--verbose``
^^^^^^^^^^^^^^^^

**Description:** Provide more verbose output if specified.

**Type:** boolean flag



``-C/--clade``
^^^^^^^^^^^^^^

**Description:** Conduct analysis on specific clade identified by CSV taxon list

**Type:** string; **Default:** None



``-L/--lnlike-thresh``
^^^^^^^^^^^^^^^^^^^^^^

**Description:** The lnlike threshhold that is the minimum value by which the log-likelihood value of the best-likelihood tree must be higher than the second-best-likelihood tree for the replicate to register as the best-likelihood topology rather than 'uncertain'. If set to zero, this turns off likelihood evaluation mode and invokes tree inference mode where a tree is simply inferred from the alignment without considering likelihood (QI values are N/A in this case).

**Type:** float; **Default:** 2.0



``-O/--min-overlap``
^^^^^^^^^^^^^^^^^^^^

**Description:** The minimum sites required to be sampled for all taxa in a given quartet.

**Type:** integer; **Default:** None



``-P/--paup``
^^^^^^^^^^^^^

**Description:** Use PAUP instead of RAxML.

**Type:** boolean flag



``-V/--verbout``
^^^^^^^^^^^^^^^^

**Description:** Provide output of the frequencies of each topology and QC.

**Type:** boolean flag



``-X/--raxml-executable``
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The name (or absolute path) of the raxml executable to be used for calculating likelihoods on quartet topologies.(default='raxml')

**Type:** None; **Default:** None


.. merge_output:

merge_output
============

Description
-----------

Combines RESULT.node.scores.csv files from separate
runs for the same phylogeny into a single set of csv and tree outputs.

http://www.github.com/FePhyFoFum/quartetsampling


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-d/--nodedata`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** file containing paths of one or moreRESULT.node.score.csv files

**Type:** None; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** new output files prefix

**Type:** None; **Default:** None



``-t/--tree`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** tree file in Newick format

**Type:** file path; **Default:** None



``-c/--clade``
^^^^^^^^^^^^^^

**Description:** ==SUPPRESS==

**Type:** None; **Default:** None



``-p/--stopk``
^^^^^^^^^^^^^^

**Description:** ==SUPPRESS==

**Type:** integer; **Default:** None



``-s/--startk``
^^^^^^^^^^^^^^^

**Description:** ==SUPPRESS==

**Type:** integer; **Default:** 0



``-v/--verbose``
^^^^^^^^^^^^^^^^

**Description:** None

**Type:** boolean flag


.. query_tree:

query_tree
==========

Description
-----------

Tree query script to find specific nodes numbers in large trees
when using the post-run annotated trees.

http://www.github.com/FePhyFoFum/quartetsampling


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-c/--clade``
^^^^^^^^^^^^^^

**Description:** ==SUPPRESS==

**Type:** None; **Default:** None



``-d/--data``
^^^^^^^^^^^^^

**Description:** CSV output from quartet_sampling (RESULT.node.score.csv)

**Type:** file path; **Default:** None



``-p/--stopk``
^^^^^^^^^^^^^^

**Description:** ==SUPPRESS==

**Type:** integer; **Default:** None



``-s/--startk``
^^^^^^^^^^^^^^^

**Description:** ==SUPPRESS==

**Type:** integer; **Default:** 0



``-t/--tree``
^^^^^^^^^^^^^

**Description:** input tree in newick format

**Type:** file path; **Default:** None



``-v/--verbose``
^^^^^^^^^^^^^^^^

**Description:** verbose screen output

**Type:** boolean flag


.. calc_qstats:

calc_qstats
===========

Description
-----------
Calculate basic statistics on the
   RESULTS.node.score.csv output file
   from quartet_sampling
   

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-d/--data`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** RESULT.node.score.csv file output fromquartet_sampling.py

**Type:** file path; **Default:** None



``-c/--clade``
^^^^^^^^^^^^^^

**Description:** specify a clade using a comma-separatedlist of 2+ descendant taxa

**Type:** None; **Default:** None



``-o/--out``
^^^^^^^^^^^^

**Description:** output file path for statistics

**Type:** file path; **Default:** None



``-p/--stopk``
^^^^^^^^^^^^^^

**Description:** stopping branch numerical index

**Type:** integer; **Default:** None



``-s/--startk``
^^^^^^^^^^^^^^^

**Description:** starting branch numerical index

**Type:** integer; **Default:** 0



``-v/--verbose``
^^^^^^^^^^^^^^^^

**Description:** verbose screen output

**Type:** boolean flag


