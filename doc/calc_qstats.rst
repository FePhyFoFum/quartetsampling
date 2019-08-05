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



``--data`` (required)
^^^^^^^^^^^^^^^^^^^^^

**Description:** RESULT.node.score.csv file output fromquartet_sampling.py

**Type:** file path; **Default:** None



``--clade``
^^^^^^^^^^^

**Description:** specify a clade using a comma-separatedlist of 2+ descendant taxa

**Type:** None; **Default:** None



``--out``
^^^^^^^^^

**Description:** output file path for statistics

**Type:** file path; **Default:** None



``--startk``
^^^^^^^^^^^^

**Description:** starting branch numerical index

**Type:** integer; **Default:** 0



``--stopk``
^^^^^^^^^^^

**Description:** stopping branch numerical index

**Type:** integer; **Default:** None



``--verbose``
^^^^^^^^^^^^^

**Description:** verbose screen output

**Type:** boolean flag


