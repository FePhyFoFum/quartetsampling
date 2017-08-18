.. _intro:

#######
Outputs
#######

RESULT.node.scores.csv
======================
A comma-separated values (CSV) file with 
* **node_label** The label of the internal branch (QS###) or terminal branch (original label)
* **freq0** The number of concordant replicates over the non-uncertain total
* **qc** The Quartet Concordance score (internal branches only; measures frequency of concordant over discordant)
* **qd** The Quartet Differential score (internal branches only; measures skew in the two discordant tree counts)
* **qi** The Quartet Informativeness score (internal branches only; measures number of replicates that fail likelihood cutoff)
* **qf** The Quartet Fidelity score (terminal branches only; measures the number of replicates for which this branch produced concordant quartets)
* **diff** Example likelihood differentinal from last replicate (used for diagnostic purposes)
* **num_replicates** Number of replicates actually sampled per branch (may not equal the number specified by ``-N`` for internal branches when fewer replicates are possible, and will not equal that number for terminal branches).
* **notes** Additional notes (enables user to add custom notes to their output after running; e.g., labeling key branches)

RESULT.labeled.tre
==================
A Newick tree with each internal branch labeled with their QS## identifier.

RESULT.labeled.tre.freq/qc/qd/qu
================================
A Newick tree with the each internal branch labeled with frequency of concordant replicates or QC/QD/QI scores.

RESULT.labeled.tre.figtree 
==========================
A FigTree format phylogeny <http://tree.bio.ed.ac.uk/software/figtree/> that contains all QS scores and a "score" field with QC/QD/QI for internal branches.

