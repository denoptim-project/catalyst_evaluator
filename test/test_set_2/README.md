# Test Set 2
This set of catalysts includes all those catalysts that were not tested experimentally, but are considered according to evaluations performed by others and published in scientific articles.

## Run the evaluations
Any evaluation starts from a graph-based definitions of the ligands. Here, we arrume a list of such candidate ligands set is available, and each such combination is saved as `*.sdf` file. These files are generated coording to [the instructions on how to make graph representations of ligand sets](../../README.md#Evaluation-of-Catalysts).

1. We first prepare the job definitions:
```
ls *sdf | while read f ; do n=$(basename $f .sdf); echo "../../evaluate_catalyst.sh -i $f" > $n.sh ; chmod +x $n.sh ; done
```

2. Then, we submit such job via a scheduler (here we use [at](https://en.wikipedia.org/wiki/At_(command))):
```
ls *.sh | while read f ; do at now -f $f ; done
```

3. Perform the analysis: __TODO__

## TODO
1. Process results of ongoing calculations to keep all the 8 molecular modes produced by each catalyst evaluation.

2. Produce statistics (i.e., update jupyter notebook)
 
