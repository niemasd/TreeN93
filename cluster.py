#!/usr/bin/env python3
'''
Given a TreeN93 tree structure, compute clusters of the leaves
'''
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
VERBOSE = False

# compute the clustering that maximizes the number of non-singleton clusters
def max_non_singleton(tree):
    best = None; num_clusters = 0; tree.resolve_polytomies()
    for node in tree.traverse_postorder(): # fix missing labels caused by resolving polytomies
        if node.is_leaf():
            node.leafdist = 0
        else:
            node.leafdist = node.children[0].leafdist + node.children[0].edge_length
    for dist,node in tree.traverse_rootdistorder(ascending=False,leaves=False):
        if node.children[0].is_leaf() and node.children[1].is_leaf(): # merging 2 singletons
            num_clusters += 1
        elif not node.children[0].is_leaf() and not node.children[1].is_leaf(): # merging 2 clusters
            num_clusters -= 1
        if best is None or num_clusters > best[1]:
            best = (float(node.leafdist),num_clusters)
    if VERBOSE:
        stderr.write("Optimal Threshold: %f\nNumber of Non-Singleton Clusters: %d\n"%best)
    clusters = list(); to_explore = Queue(); to_explore.put(tree.root)
    while not to_explore.empty():
        node = to_explore.get()
        if node.is_leaf():
            clusters.append([str(node)])
        elif float(node.leafdist) > best[0]:
            for c in node.children:
                to_explore.put(c)
        else:
            clusters.append([str(u) for u in node.traverse_leaves()])
    return clusters

MODES = {
    'max_non_singleton': max_non_singleton
}

# main execution
if __name__ == "__main__":
    from treeswift import read_tree_newick
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    parser.add_argument('-m', '--mode', required=False, type=str, default='max_non_singleton', help="Clustering Mode")
    parser.add_argument('-v', '--verbose', action="store_true", help="Verbose Mode")
    args = parser.parse_args()
    if args.mode.lower() not in MODES:
        raise ValueError("Invalid mode: %s. Valid options: %s" % (args.mode, ', '.join(MODES.keys())))
    if args.verbose:
        VERBOSE = True; global stderr; from sys import stderr; from time import time
    if VERBOSE:
        stderr.write("=== INPUTS ===\n")
        stderr.write("Input File:  %s\n"%args.input)
        stderr.write("Output File: %s\n"%args.output)
        stderr.write("Mode:        %s\n"%args.mode)
        stderr.write('\n'); stderr.flush()
    if args.input == 'stdin':
        from sys import stdin; treestr = stdin.read().strip()
    elif args.input.lower().endswith('.gz'):
        from gzip import open as gopen; treestr = gopen(args.input).read().decode().strip()
    else:
        treestr = open(args.input).read().strip()
    if args.output == 'stdout':
        from sys import stdout as outfile
    elif args.output.lower().endswith('.gz'):
        from gzip import open as gopen; outfile = gopen(args.output,'w')
    else:
        outfile = open(args.output,'w')

    # compute and output clusters
    if VERBOSE:
        stderr.write("=== RESULTS ===\n")
        stderr.write("Computing transmission clusters...\n"); START = time()
    clusters = MODES[args.mode](read_tree_newick(treestr))
    cluster_num = 1; outfile.write("SequenceName\tClusterNumber\n")
    for c in clusters:
        if len(c) == 1:
            outfile.write(c[0]); outfile.write('\t-1\n'); continue
        for u in c:
            outfile.write(u); outfile.write('\t%d\n'%cluster_num)
        cluster_num += 1
    outfile.close()
    if VERBOSE:
        END = time(); stderr.write("Total runtime: %s seconds\n" % (END-START))
