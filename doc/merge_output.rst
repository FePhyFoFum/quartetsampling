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


