.. _output:

#######
Outputs
#######

RESULT.node.scores.csv
======================
A comma-separated values (CSV) file with: 

* **node_label** The label of the internal branch (QS###) or terminal branch (original label)
* **freq0** The number of concordant replicates over the non-uncertain total
* **qc** The Quartet Concordance score (internal branches only; measures frequency of concordant over discordant)
* **qd** The Quartet Differential score (internal branches only; measures skew in the two discordant tree counts)
* **qi** The Quartet Informativeness score (internal branches only; measures number of replicates that fail likelihood cutoff)
* **qf** The Quartet Fidelity score (terminal branches only; measures the number of replicates for which this branch produced concordant quartets)
* **diff** Example likelihood differential from last replicate (used for diagnostic purposes)
* **num_replicates** Number of replicates actually sampled per branch (may not equal the number specified by ``-N`` for internal branches when fewer replicates are possible, and will not equal that number for terminal branches).
* **notes** Additional notes (enables user to add custom notes to their output after running; e.g., labeling key branches)

RESULT.node.counts.csv
======================
A tab-separated values (TSV) file with:

* **node_label** The label of the internal branch (QS###) or terminal branch (original label)
* **count0** The count of the number of QS replicates for the concordant quartet arrangement.
* **count1** The count of the number of QS replicates for one of the discordant quartet arrangements.
* **count2** The count of the number of QS replicates for the other discordant quartet arrangement.
* **topo0**, **topo1**, and **topo2** provide example quartet trees (in parenthetical notation) showing the arragement of a representative quartet of taxa spanning the node for the arragements corresponding to the counts.

This file allows you to check at a node which discordant option is more common in cases where QD is low, indicating higher presence of one discordant option.

RESULT.labeled.tre
==================
A Newick tree with each internal branch labeled with their QS## identifier.

RESULT.labeled.tre.freq/qc/qd/qu
================================
A Newick tree with the each internal branch labeled with frequency of concordant replicates or QC/QD/QI scores.

RESULT.labeled.tre.figtree 
==========================
A FigTree format phylogeny <http://tree.bio.ed.ac.uk/software/figtree/> that contains all QS scores and a "score" field with QC/QD/QI for internal branches.

