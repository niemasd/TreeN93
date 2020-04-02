#!/usr/bin/env python3
'''
TreeN93: Non-parametric transmission clustering from pairwise phylogenetic distances (Niema Moshiri 2018)
'''
from treeswift import Node
from niemads import DisjointSet
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
VERBOSE = False

# parse TN93 input
def parse_tn93(infile):
    dists = list() # list of (d,u,v) tuples
    for line in infile:
        if isinstance(line,bytes):
            u,v,d = line.decode().strip().split(',')
        else:
            u,v,d = line.strip().split(',')
        if d == 'Distance':
            continue
        dists.append((float(d),u,v))
    dists.sort(); return dists

# create TreeN93 tree from input pairwise distances
def dist_to_tree(dists,missing):
    dists.sort()
    if missing < dists[-1][0]:
        raise ValueError("Missing value must be larger than max distance")
    root = dict() # root[s] = the root of the tree corresponding to set s
    ds = DisjointSet()
    for d,u,v in dists:
        if u not in ds:
            ds.add(u); root[u] = Node(label=u)
        if v not in ds:
            ds.add(v); root[v] = Node(label=v)
        su = ds.find(u); sv = ds.find(v)
        if su != sv:
            ru = root[su]; del root[su]; rv = root[sv]; del root[sv]; ds.union(su,sv)
            p = Node(label=d); p.add_child(ru); p.add_child(rv); root[ds.find(su)] = p
            if ru.is_leaf():
                ru.edge_length = d
            else:
                ru.edge_length = d - ru.label
            if rv.is_leaf():
                rv.edge_length = d
            else:
                rv.edge_length = d - rv.label
    if VERBOSE:
        stderr.write("Number of Individuals: %d\n"%len(ds))
    if len(root.values()) > 1:
        toproot = Node(label=missing)
        for r in root.values():
            toproot.add_child(r); r.edge_length = missing - r.label
        return [toproot]
    return list(root.values())[0]

# main execution
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    parser.add_argument('-m', '--missing', required=False, type=float, default=float('inf'), help="Value for Missing Distances")
    parser.add_argument('-v', '--verbose', action="store_true", help="Verbose Mode")
    args = parser.parse_args()
    if args.verbose:
        VERBOSE = True; global stderr; from sys import stderr; from time import time
    if VERBOSE:
        stderr.write("=== INPUTS ===\n")
        stderr.write("Input File:  %s\n"%args.input)
        stderr.write("Output File: %s\n"%args.output)
        stderr.write('\n'); stderr.flush()
    if args.input == 'stdin':
        from sys import stdin; infile = stdin.read().strip().splitlines()
    elif args.input.lower().endswith('.gz'):
        from gzip import open as gopen; infile = gopen(args.input).read().decode().strip().splitlines()
    else:
        infile = open(args.input).read().strip().splitlines()

    # compute and output TreeN93 tree
    if VERBOSE:
        stderr.write("=== COMPUTE TREEN93 TREE ===\n"); START = time()
    dists = parse_tn93(infile)
    root = dist_to_tree(dists,args.missing)
    if args.output == 'stdout':
        print(root.newick())
    elif args.output.lower().endswith('.gz'):
        f = gopen(args.output,'wb',9); f.write((('%s;\n' % root.newick()).encode())); f.close()
    else:
        f = open(args.output,'w'); f.write('%s;\n' % root.newick()); f.close()
    if VERBOSE:
        END = time(); stderr.write("Total runtime: %s seconds\n" % (END-START))
