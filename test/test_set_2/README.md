# Test Set 2
This set of catalysts includes all those catalysts that were not tested experimentally, but are considered according to evaluations performed by others and published in scientific articles.

## Run Catalysts Evaluations
Any evaluation starts from a graph-based definitions of the ligands. Here, we assume a list of such candidate ligand sets is available, and each such combination is saved as `*.sdf` file. These files are generated according to [the instructions on how to make graph representations of ligand sets](../../README.md#Evaluation-of-Catalysts).

1. We first prepare the job definitions:
```
ls *sdf | while read f ; do n=$(basename $f .sdf); echo "../../evaluate_catalyst.sh -i $f" > $n.sh ; chmod +x $n.sh ; done
```

2. Then, we submit such job via a scheduler (here we use [at](https://en.wikipedia.org/wiki/At_(command))):
```
ls *.sh | while read f ; do at now -f $f ; done
```

3. Remove tmp files (WARNING: you may want to take a backup copy!)

```
origin=$(pwd); ls */*_out.sdf | while read f ; do d=$(dirname $f); if ! grep -q FITNESS $f ; then echo "ERROR: missing fitness for $d" ; continue ; fi; dd=${d}_tmp ; mkdir $dd ; cp $d/${d}_out.sdf $d/${d}.tar.gz $d/${d}.sdf  $dd ; mkdir -p trash/${d}_trash ; cd $d ; cp -r * ../trash/${d}_trash ; rm -r * ; mv ../$dd/* . ; tar -xzvf $d.tar.gz $d/${d}_outSubPreDFT-D.sdf $d/${d}_outSubPreDFT-X.sdf $d/${d}_outSubPreDFT-Z.sdf $d/${d}_outSubXTB-A.sdf $d/${d}_outSubXTB-C.sdf $d/${d}_outSubXTB-E.sdf $d/${d}_outSubXTB-F.sdf $d/${d}_outSubXTB-L.sdf ; mv $d/* . ; rm -r "$d" ${d}.tar.gz ; rm -r ../$dd ; cd $origin ; done
```

4. Perform the analysis. Run the jupyter notebook [Analysis_test_set_2.ipynb](Analysis_test_set_2.ipynb).


## File Naming
The files are named according to the following internal convention, which is encoded in the scripts provided under [../../src/](../../src/).

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
