#!/usr/bin/env python3
'''
TreeN93: Non-parametric transmission clustering from pairwise phylogenetic distances (Niema Moshiri 2018)
'''
from treeswift import Node
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
VERBOSE = False

# helper disjoint set class
class DisjointSet:
    def __init__(self): # initialize
        self.parent = dict() # parent[u] = parent of node u
        self.num_below = dict() # num_below[u] = number of nodes below u (including u) (only current for sentinels)
    def __contains__(self,x):
        return x in self.parent
    def __len__(self): # number of elements in Disjoint Set
        return len(self.parent)
    def add(self,x): # add x as a sentinel node
        if x in self:
            raise ValueError("Node already exists: %s"%x)
        self.parent[x] = None; self.num_below[x] = 1
    def remove(self,x): # remove x from Disjoint Set
        if x not in self:
            raise ValueError("Node not found: %s"%x)
        children = [n for n in self.parent if self.parent[n] == x]
        if len(children) != 0:
            if self.parent[x] is None:
                p = children.pop(); self.parent[p] = None; self.num_below[p] = self.num_below[x] - 1
            else:
                p = self.parent[x]; self.num_below[p] -= 1
            for c in children:
                self.parent[c] = p
        del self.parent[x]; del self.num_below[x]
    def find(self,x): # return the sentinel node of x
        if x not in self:
            raise ValueError("Node not found: %s"%x)
        explored = Queue(); curr = x
        while self.parent[curr] is not None:
            explored.put(curr); curr = self.parent[curr]
        while not explored.empty():
            self.parent[explored.get()] = curr # path compression
        return curr
    def union(self,x,y): # union the sets containing x and y
        if x not in self:
            raise ValueError("Node not found: %s"%x)
        if y not in self:
            raise ValueError("Node not found: %s"%y)
        sx = self.find(x); sy = self.find(y)
        if self.num_below[sx] > self.num_below[sy]:
            self.parent[sy] = sx; self.num_below[sx] += self.num_below[sy]
        else:
            self.parent[sx] = y; self.num_below[sy] += self.num_below[sx]

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
    if args.output == 'stdout':
        from sys import stdout as outfile
    elif args.output.lower().endswith('.gz'):
        from gzip import open as gopen; outfile = gopen(args.output,'w')
    else:
        outfile = open(args.output,'w')

    # compute and output TreeN93 tree
    if VERBOSE:
        stderr.write("=== COMPUTE TREEN93 TREE ===\n"); START = time()
    dists = parse_tn93(infile)
    root = dist_to_tree(dists,args.missing)
    outfile.write(root.newick()); outfile.write(';\n'); outfile.close()
    if VERBOSE:
        END = time(); stderr.write("Total runtime: %s seconds\n" % (END-START))
