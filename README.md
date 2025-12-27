(A) Parsimony analyses are implemented in the file "leptop.tnt". This will run analyses under implied weighting, using three concavities. Simply delete the "p/" in line 127.

(B) Non-clock Bayesian analyses are called "leptop_noclockstate.nex" and "leptop_noclocksource.nex"

(C) Clock Bayesian analyses are in the folder "clock". Files are named such that they inform the clock, sampling strategy, and partitioning type:

- clock: "igr" or "tk02"
- sampling: "fossiltip", "ran" (random), "div" (diversity)
- partitioning: "state" (by the number of states), "source" (vegetative vs fertile)

For instance, the file "leptop_igrdivsource.nex" implemented an scheme with igr clock, diversity sampling, and by-source partitioning. 

(D) In Zenodo (10.5281/zenodo.18064595), you will find an R script to replicate analyses in Urrea et al. (2026) and its input files. It includes code for:

- Computing Disparity through time, using dispRity (Guillerme, 2018) and Claddis (Lloyd, 2016) - Step 8
- Estimating natural selection modes and strength (EvoPhylo, Simões et al., 2023), and processing and plotting output from MrBayes (ggtree) - Steps 10-11
