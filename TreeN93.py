#!/usr/bin/env python3
'''
TreeN93: Non-parametric transmission clustering from pairwise phylogenetic distances (Niema Moshiri 2018)
'''
from treeswift import Node,read_tree_newick
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

# create TreeN93 tree from input TN93 distance file
def dist_to_tree(dists,missing):
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
    return list(root.values())

# compute the clustering that maximizes the number of non-singleton clusters
def tree_to_clusters(root):
    best = None; num_clusters = 0
    for dist,node in root.traverse_rootdistorder(ascending=False,leaves=False):
        assert len(node.children) == 2, "TreeN93 tree not fully bifurcating"
        if node.children[0].is_leaf() and node.children[1].is_leaf(): # merging 2 singletons
            num_clusters += 1
        elif not node.children[0].is_leaf() and not node.children[1].is_leaf(): # merging 2 clusters
            num_clusters -= 1
        if best is None or num_clusters > best[1]:
            best = (float(node.label),num_clusters)
    if VERBOSE:
        stderr.write("Optimal Threshold: %f\nNumber of Non-Singleton Clusters: %d\n"%best)
    clusters = list(); to_explore = Queue(); to_explore.put(root)
    while not to_explore.empty():
        node = to_explore.get()
        if node.is_leaf():
            clusters.append([str(node)])
        elif float(node.label) > best[0]:
            for c in node.children:
                to_explore.put(c)
        else:
            clusters.append([str(u) for u in node.traverse_leaves()])
    return clusters

# main execution
if __name__ == "__main__":
    from subprocess import run
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File")
    parser.add_argument('-t', '--num_threads', required=False, type=int, default=1, help="Number of Threads")
    parser.add_argument('-m', '--missing', required=False, type=float, default=float('inf'), help="Value for Missing Distances")
    parser.add_argument('-mafft', '--mafft_path', required=False, type=str, default='mafft', help="MAFFT Executable Path")
    parser.add_argument('-iqtree', '--iqtree_path', required=False, type=str, default='iqtree', help="IQ-TREE Executable Path")
    parser.add_argument('-iqm', '--iqtree_model', required=False, type=str, default='MFP', help="IQ-Tree Model")
    parser.add_argument('-v', '--verbose', action="store_true", help="Verbose Mode")
    args = parser.parse_args()
    if args.verbose:
        VERBOSE = True; global stderr; from sys import stderr; from time import time
    if VERBOSE:
        stderr.write("=== INPUTS ===\n")
        stderr.write("Input File:     %s\n"%args.input)
        stderr.write("# Threads:      %d\n"%args.num_threads)
        stderr.write("MAFFT Path:     %s\n"%args.mafft_path)
        stderr.write("IQ-TREE Path:   %s\n"%args.iqtree_path)
        stderr.write("IQ-TREE Model:  %s\n"%args.iqtree_model)
        stderr.write('\n'); stderr.flush()
    if args.input == 'stdin':
        from sys import stdin as infile
    elif args.input.lower().endswith('.gz'):
        from gzip import open as gopen; infile = gopen(args.input)
    else:
        infile = open(args.input)

    # compute and output TreeN93 trees
    if VERBOSE:
        stderr.write("=== COMPUTE TREEN93 TREES ===\n"); START = time()
    infile_lines = infile.read()
    if isinstance(infile_lines, bytes):
        infile_lines = infile_lines.decode().strip().splitlines()
    else:
        infile_lines = infile_lines.strip().splitlines()
    if infile_lines[0][-1] == ';': # TreeN93 Tree (Newick)
        if VERBOSE:
            stderr.write("Input file parsed as TreeN93 tree(s) (Newick)\n")
        tree_roots = [read_tree_newick(l).root for l in infile_lines]
    else:
        if infile_lines[0][0] == '>': # Sequences (FASTA)
            if VERBOSE:
                stderr.write("Input file parsed as sequences (FASTA)\n")
            from tempfile import NamedTemporaryFile
            tmp = NamedTemporaryFile(mode='w')
            for l in infile_lines:
                tmp.write(l); tmp.write('\n')
            tmp.flush()
            if VERBOSE:
                stderr.write("Running MAFFT...\n"); stderr.flush()
            o = run([args.mafft_path, '--auto', '--thread', str(args.num_threads), tmp.name], capture_output=True)
            f = open('%s.mafft.aln'%args.input,'wb'); f.write(o.stdout); f.close()
            f = open('%s.mafft.log'%args.input,'wb'); f.write(o.stderr); f.close()
            if VERBOSE:
                stderr.write("Running IQ-TREE...\n"); stderr.flush()
            run([args.iqtree_path, '-s', '%s.mafft.aln'%args.input, '-pre', '%s.iqtree'%args.input, '-nt', str(args.num_threads), '-ntmax', str(args.num_threads), '-redo'], capture_output=True)
            if VERBOSE:
                stderr.write("Computing pairwise distances...\n"); stderr.flush()
            dm = read_tree_newick('%s.iqtree.treefile'%args.input).distance_matrix()
            nodes = list(dm.keys())
            dists = [(dm[nodes[i]][nodes[j]], nodes[i].label, nodes[j].label) for i in range(len(nodes)-1) for j in range(i+1,len(nodes))]
            dists.sort()
        else: # Pairwise Distances
            if VERBOSE:
                stderr.write("Input file parsed as pairwise distances\n")
            dists = parse_tn93(infile_lines)
        tree_roots = dist_to_tree(dists,args.missing)
        outfile_trees = open('%s.treen93.nwk'%args.input,'w')
        for root in tree_roots:
            outfile_trees.write(root.newick()); outfile_trees.write(';\n')
        outfile_trees.close()
        if VERBOSE:
            stderr.write('\n'); stderr.flush()

    # compute and output clusters maximizing number of non-singleton clusters
    if VERBOSE:
        stderr.write("=== RESULTS ===\n")
        stderr.write("Computing transmission clusters...\n")
    clusters = list()
    for root in tree_roots:
        for cluster in tree_to_clusters(root):
            clusters.append(cluster)
    cluster_num = 1
    outfile_clusters = open('%s.clusters.tsv'%args.input,'w')
    outfile_clusters.write("SequenceName\tClusterNumber\n")
    for c in clusters:
        if len(c) == 1:
            outfile_clusters.write(c[0]); outfile_clusters.write('\t-1\n'); continue
        for u in c:
            outfile_clusters.write(u); outfile_clusters.write('\t%d\n'%cluster_num)
        cluster_num += 1
    outfile_clusters.close()
    if VERBOSE:
        END = time(); stderr.write("Total runtime: %s seconds" % (END-START))
