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



``--align/--alignment`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Alignment file in "relaxed phylip" format, as used by RAxML.

**Type:** file path; **Default:** None



``--reps/--number-of-reps`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The number of replicate quartet topology searches to be performed at each node.

**Type:** integer; **Default:** 100



``--threads/--number-of-threads`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The number of parallel threads to be used by Python for quartet topology searches.

**Type:** integer; **Default:** 1



``--tree`` (required)
^^^^^^^^^^^^^^^^^^^^^

**Description:** The input tree in Newick (parenthetical) format.

**Type:** file path; **Default:** None



``--calc-qdstats``
^^^^^^^^^^^^^^^^^^

**Description:** EXPERIMENTAL: Calculates Chi-square test for QD tree frequencies. Use only  if Scipy is available. Will increase running time.

**Type:** boolean flag



``--clade``
^^^^^^^^^^^

**Description:** Conduct analysis on specific clade identified by CSV taxon list

**Type:** string; **Default:** None



``--data-type``
^^^^^^^^^^^^^^^

**Description:** (nuc)leotide, (amino) acid, or (cat)egorical data

**Type:** None; **Default:** ['nuc']

**Choices:** ('nuc', 'amino', 'cat')


``--engine``
^^^^^^^^^^^^

**Description:** Name of the program to use to infer trees or evaluate tree model likelihoods.

**Type:** None; **Default:** ('raxml-ng',)

**Choices:** ('raxml-ng', 'raxml', 'paup', 'iqtree')


``--engine-exec``
^^^^^^^^^^^^^^^^^

**Description:** Full file path of the tree inference or likelihood evaluation engine.

**Type:** None; **Default:** None



``--engine-model``
^^^^^^^^^^^^^^^^^^

**Description:** Advanced: specify a custom model name for the tree engine

**Type:** None; **Default:** None



``--genetrees``
^^^^^^^^^^^^^^^

**Description:** Use partitions file (RAxML format) to divide the alignment into separate gene tree regions. Gene alignments will be sampled random for the quartet topology searches.

**Type:** file path; **Default:** None



``--ignore-errors``
^^^^^^^^^^^^^^^^^^^

**Description:** Ignore RAxML and PAUP erroneous runs

**Type:** boolean flag



``--lnlike/--lnlike-thresh``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The lnlike threshhold that is the minimum value by which the log-likelihood value of the best-likelihood tree must be higher than the second-best-likelihood tree for the replicate to register as the best-likelihood topology rather than 'uncertain'. If set to zero, this turns off likelihood evaluation mode and invokes tree inference mode where a tree is simply inferred from the alignment without considering likelihood (QI values are N/A in this case).

**Type:** float; **Default:** 2.0



``--low-mem``
^^^^^^^^^^^^^

**Description:** Do not store large alignment in memory for whole-alignment (non-genetree) mode

**Type:** boolean flag



``--max-random-sample-proportion``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** The proportion of possible replicates explored unsuccessfully by the random generation procedure before it gives up. Because this generates random replicates, it takes progressively longer as it proceeds. To avoid long runtimes, the recommended range is < 0.5 (which is the default).

**Type:** float; **Default:** None



``--min-overlap``
^^^^^^^^^^^^^^^^^

**Description:** The minimum sites required to be sampled for all taxa in a given quartet.

**Type:** integer; **Default:** None



``--partitions``
^^^^^^^^^^^^^^^^

**Description:** Partitions file in RAxML format. If omitted then the entire alignment will be treated as one partition for all quartet replicate topology searches.

**Type:** file path; **Default:** None



``--result-prefix``
^^^^^^^^^^^^^^^^^^^

**Description:** A prefix to put on the result files.

**Type:** string; **Default:** None



``--results-dir``
^^^^^^^^^^^^^^^^^

**Description:** A directory to which output files will be saved. If not supplied, the current working directory will be used. (default is current folder).

**Type:** file path; **Default:** None



``--retain-temp``
^^^^^^^^^^^^^^^^^

**Description:** Do not remove temporary files

**Type:** boolean flag



``--start-node-number``
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** An integer denoting the node to which to start from. Nodes will be read from topologically identical (and isomorphic!) input trees in deterministic order, so this argument may be  used to restart at an intermediate position (in case the previous run was canceled before completion, for example).

**Type:** integer; **Default:** None



``--stop-node-number``
^^^^^^^^^^^^^^^^^^^^^^

**Description:** An integer denoting the node at which to stop. Will include nodes with indices <= the stop node number. This argument may be used to limit the length of a given run in case only a certain part of the tree is of interest. Nodes will be read from topologically identical (and isomorphic!) input trees in deterministic order.

**Type:** integer; **Default:** None



``--temp-dir``
^^^^^^^^^^^^^^

**Description:** A directory to which temporary files will be saved. If not supplied, 'QuartetSampling' will be created in the current working directory. When specifying a custom temporary output the characters 'QuartetSampling' must appear in the directory name to prevent accidental file deletion. (default='./QuartetSampling'

**Type:** file path; **Default:** None



``--verbose``
^^^^^^^^^^^^^

**Description:** Provide more verbose output if specified.

**Type:** boolean flag



``--verbout``
^^^^^^^^^^^^^

**Description:** Provide output of the frequencies of each topology and QC.

**Type:** boolean flag


