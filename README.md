# TreeN93
TreeN93 is a non-parametric Python 3 tool that performs transmission cluster identification from pairwise TN93 distances.

## Clustering Objective

Let *C*(*S*|*t*) denote a clustering of the sequences in a set *S* using a distance threshold *t* such that, for all pairs of elements (*u*,*v*) in *S*, if the TN93 distance between *u* and *v* is below *t*, *u* and *v* are placed in the same cluster, and the number of clusters is maximized. Let *N*(*S*|*t*) denote the number of clusters containing more than one element in *C*(*S*|*t*). TreeN93 outputs a clustering *C*(*S*|*t'*) that maximizes *N*(*S*|*t'*) over all possible *t'*. In other words, if you were to run HIV-TRACE using all possible thresholds *t* and then pick the threshold *t'* that yielded the clustering with the most non-singleton clusters, that would be the clustering output by TreeN93.

## Installation
TreeN93 is written in Python 3 and depends on the [TreeSwift](https://github.com/niemasd/TreeSwift) Python package. Once TreeSwift is installed, simply download [TreeN93.py](https://github.com/niemasd/TreeN93/blob/master/TreeN93.py) to your machine and make it executable.

## Usage
TreeN93 can be used as follows:

```bash
usage: TreeN93.py [-h] [-i INPUT] -o OUTPRE [-v]

optional arguments:
  -h, --help                    show this help message and exit
  -i INPUT, --input INPUT       Input TN93 File (default: stdin)
  -o OUTPRE, --outpre OUTPRE    Output Prefix (default: None)
  -v, --verbose                 Print Verbose Messages to Standard Error (default: False)
```

The input file must be in the format output by [tn93](https://github.com/veg/tn93). Note that TreeN93 requires *all* pairwise distances to be output, whereas tn93 only outputs distances below the chosen threshold (0.015, or 1.5%, by default) and with an overlap longer than the chosen threshold (100 bases by default), so you will want to run it using `-t 1` to specify a distance threshold of 1 (i.e., 100%) and `-l 1` to specify an overlap threshold of 1. For example:

```bash
tn93 -t 1 -l 1 my_sequences.fas > my_sequences.tn93.csv
TreeN93.py -i my_sequences.tn93.csv -o my_sequences
```

Note that TreeN93 can also read TN93 distances from standard input, so it can be easily piped into an existing workflow. For example:

```bash
zcat my_sequences.fas.gz | tn93 -t 1 -l 1 | TreeN93.py -o my_sequences
```

## Output Files
TreeN93 outputs two files: a "clustering" file (`prefix.clusters.txt`) and a TreeN93 tree structure (`prefix.trees.nwk`). The "clustering" file is the resulting clustering *C*(*S*|*t'*) as described above, output in the same format as [TreeCluster](https://github.com/niemasd/TreeCluster). The TreeN93 tree structure is a Newick-format tree that, if cut *t* distance above the leaves, will yield the HIV-TRACE clustering obtained using a distance threshold of *t*. For example, say you want to obtain TreeN93 clusters but you would also like to obtain the clusters that HIV-TRACE would have found with its default threshold of 0.015, you could do so as follows:

```bash
zcat my_sequences.fas.gz | tn93 -t 1 -l 1 | TreeN93.py -o my_sequences
TreeCluster.py -i my_sequences.trees.nwk -m leaf_dist_min -t 0.015 -o my_sequences.hivtrace.txt.gz
```

If the input TN93 distance file does not actually contain all *n*(*n*-1)/2 pairwise distances, the tree structure file will contain one tree for each component of the graph *G* = (*V*,*E*) where *V* contains a vertex *v* for each sequence and *E* contains an edge (*u*,*v*,*d*) for each pairwise distance *d* between sequences *u* and *v* in the input, and the clustering file will contain the results of running TreeN93 on each component.
