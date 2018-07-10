#!/usr/bin/env python3
'''
TreeN93: Construct a hierarchical tree from TN93 distances
'''
from treeswift import Node
try:
    from Queue import Queue
except ImportError:
    from queue import Queue

# helper disjoint set class
class DisjointSet:
    try:
        from Queue import Queue
    except ImportError:
        from queue import Queue
    def __init__(self): # initialize
        self.parent = dict() # parent[u] = parent of node u
        self.num_below = dict() # num_below[u] = number of nodes below u (including u) (only current for sentinels)
    def __contains__(self,x):
        return x in self.parent
    def add(self,x): # add x as a sentinel node
        if x in self:
            raise ValueError("Node already exists: %s"%x)
        self.parent[x] = None; self.num_below[x] = 1
    def find(self,x): # return the sentinel node of x
        if x not in self:
            raise ValueError("Node not found: %s"%x)
        explored = Queue(); curr = x
        while self.parent[curr] != None:
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

# create TreeN93 tree from input TN93 distance file
def dist_to_tree(dists):
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
    return ['%s;'%r.newick() for r in root.values()]

# main execution
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input TN93 File")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Newick File")
    args = parser.parse_args()
    if args.input == 'stdin':
        from sys import stdin as infile
    elif args.input.lower().endswith('.gz'):
        from gzip import open as gopen; infile = gopen(args.input)
    else:
        infile = open(args.input)
    if args.output == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output,'w')
    for tree in dist_to_tree(parse_tn93(infile)):
        outfile.write(tree); outfile.write('\n')
    outfile.close()
