# TreeN93
TreeN93 is a non-parametric Python 3 tool that performs transmission cluster identification from pairwise TN93 distances.

## Installation
TreeN93 is written in Python 3 depends on the [TreeSwift](https://github.com/niemasd/TreeSwift) Python package. Once TreeSwift is installed, simply download [TreeN93.py](https://github.com/niemasd/TreeN93/blob/master/TreeN93.py) to your machine and make it executable.

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
tn93 -t 1 -l 1 my_sequences.fas | TreeN93.py -o my_sequences
```
