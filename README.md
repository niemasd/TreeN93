# TreeN93
TreeN93 is a non-parametric Python 3 tool that performs transmission cluster identification from pairwise phylogenetic (e.g. TN93) distances.

## Clustering Objective
Let *C*(*S*|*t*) denote a clustering of the sequences in a set *S* using a distance threshold *t* such that, for all pairs of elements (*u*,*v*) in *S*, if the phylogenetic distance between *u* and *v* is below *t*, *u* and *v* are placed in the same cluster, and the number of clusters is maximized. Let *N*(*S*|*t*) denote the number of clusters containing more than one element in *C*(*S*|*t*). TreeN93 outputs a clustering *C*(*S*|*t'*) that maximizes *N*(*S*|*t'*) over all possible *t'*. In other words, if you were to run HIV-TRACE using all possible thresholds *t* and then pick the threshold *t'* that yielded the clustering with the most non-singleton clusters, that would be the clustering output by TreeN93.

## Web App
TreeN93 has recently been developed into a web app that runs fully client-side in the user's web browser (i.e., the user's data are not sent anywhere) and that is much more user-friendly. We highly recommend utilizing the web app if your dataset is small enough to run fast enough, and we only recommend utilizing the command-line Python tool if your dataset is too large to run reasonably fast in the web browser.

https://github.com/Niema-Lab/TreeN93

## Installation
TreeN93 is written in Python 3 and depends on the [TreeSwift](https://github.com/niemasd/TreeSwift) and [NiemaDS](https://github.com/niemasd/NiemaDS) Python packages. Once they are installed, simply download [TreeN93.py](https://github.com/niemasd/TreeN93/blob/master/TreeN93.py) to your machine and make it executable.

## Usage
[TreeN93.py](TreeN93.py) can be used as follows:

```bash
usage: TreeN93.py [-h] [-i INPUT] [-o OUTPUT] [-m MISSING] [-v]

optional arguments:
  -h, --help                      show this help message and exit
  -i INPUT, --input INPUT         Input File (default: stdin)
  -o OUTPUT, --output OUTPUT      Output File (default: stdout)
  -m MISSING, --missing MISSING   Value for Missing Distances (default: inf)
  -v, --verbose                   Verbose Mode (default: False)
```

[TreeN93_cluster.py](TreeN93_cluster.py) can be used as follows:

```bash
usage: TreeN93_cluster.py [-h] [-i INPUT] [-o OUTPUT] [-m MODE] [-v]

optional arguments:
  -h, --help                      show this help message and exit
  -i INPUT, --input INPUT         Input File (default: stdin)
  -o OUTPUT, --output OUTPUT      Output File (default: stdout)
  -m MODE, --mode MODE            Clustering Mode (default: max_non_singleton)
  -v, --verbose                   Verbose Mode (default: False)
```

## Output Files
[TreeN93.py](TreeN93.py) outputs a TreeN93 tree structure, which is a Newick-format tree that, if cut *t* distance above the leaves, will yield the HIV-TRACE clustering obtained using a distance threshold of *t*. For example, say you want to obtain TreeN93 clusters but you would also like to obtain the clusters that HIV-TRACE would have found with its default threshold of 0.015, you could do so as follows:

```bash
zcat my_sequences.fas.gz | tn93 -t 1 -l 1 | TreeN93.py -o my_sequences.treen93.nwk
TreeCluster.py -i my_sequences.treen93.nwk -m leaf_dist_min -t 0.015
```

If the input distance file does not actually contain all *n*(*n*-1)/2 pairwise distances, the tree structure file will contain one tree for each component of the graph *G* = (*V*,*E*) where *V* contains a vertex *v* for each sequence and *E* contains an edge (*u*,*v*,*d*) for each pairwise distance *d* between sequences *u* and *v* in the input, and the clustering file will contain the results of running TreeN93 on each component.

[TreeN93_cluster.py](TreeN93_cluster.py) outputs the resulting clustering *C*(*S*|*t'*) as described above, output in the same format as [TreeCluster](https://github.com/niemasd/TreeCluster).

## Input: Pairwise Distances
Pairwise distances must be given in the same format as [tn93](https://github.com/veg/tn93). Note that TreeN93 requires *all* pairwise distances to be output, whereas tn93 only outputs distances below the chosen threshold (0.015, or 1.5%, by default) and with an overlap longer than the chosen threshold (100 bases by default), so you will want to run it using `-t 1` to specify a distance threshold of 1 (i.e., 100%) and `-l 1` to specify an overlap threshold of 1. For example:

```bash
tn93 -t 1 -l 1 my_sequences.fas > my_sequences.tn93.csv
TreeN93.py -i my_sequences.tn93.csv
```

Note that TreeN93 can also read TN93 distances from standard input, so it can be easily piped into an existing workflow. For example:

```bash
zcat my_sequences.fas.gz | tn93 -t 1 -l 1 | TreeN93.py
```

## Acknowledgements
Much thanks to [Sergei Pond](http://spond.github.io/CV.js/cv.html), [Steven Weaver](http://www.stevenweaver.org/), and [Joel Wertheim](http://id.ucsd.edu/faculty/wertheim.shtml) for their excellent work on HIV-TRACE (namely the [tn93](https://github.com/veg/tn93) component).
