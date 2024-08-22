# Test Set 1
This set of catalysts includes all those catalysts that were synthesised and tested experimentally by our labs.
The figures published in the paper can be reproduced using the Jupyter notebook [Analysis_MS_figures.ipynb](Analysis_MS_figures.ipynb).
A more general-purpose analysis of the AF scores can be performed with [Analysis_test_set_1.ipynb](Analysis_test_set_1.ipynb).


## File Naming
The files are named according to the following internal convention, which is encoded in the scripts provided under [../../src/](../../src/). Here, `name` is a string that identifies the input to the [catalyst evaluator](../../evaluate_catalyst.sh).

* `name.sdf` contains the graph-based definition of the ligand set; This is the very input to the [catalyst evaluator](../../evaluate_catalyst.sh).
* `name_out.sdf` is the final output of the catalyst evaluator; it contains the values of the fitness score, descriptors, and weights.
* `name_outSubPreDFT-D.sdf` contains the geometry of the metallacyclobutane intermediate (DFT-based optimization).
* `name_outSubPreDFT-X.sdf` contains the transition state model for productive metathesis reaction (DFT-based constrained optimization).
* `name_outSubPreDFT-Z.sdf` contains the transition state model for beta-H elimination reaction (DFT-based constrained optimization).
* `name_outSubXTB-A.sdf` contains the geometry of the catalyst precursor (xTB-based optimization).
* `name_outSubXTB-C.sdf` contains the geometry of the metallacyclobutane intermediate (xTB-based optimization).
* `name_outSubXTB-E.sdf` (Optional: produced only when X ligands differ from Cl) contains the geometry of the catalyst precursor with trans-X configuration and X=Cl (xTB-based optimization).
* `name_outSubXTB-F.sdf` contains the geometry of the catalyst precursor with cis-X configuration (xTB-based optimization).
* `name_outSubXTB-L.sdf` contains the geometry of the free dative ligand (xTB-based optimization).
