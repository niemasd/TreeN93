#!/usr/bin/env python3
'''
Compute pairwise distances using the MAFFT + IQ-TREE workflow
'''
VERBOSE = False

# main execution
if __name__ == "__main__":
    from subprocess import PIPE,run
    from treeswift import read_tree_newick
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    parser.add_argument('-t', '--num_threads', required=False, type=int, default=1, help="Number of Threads")
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
        stderr.write("Output File:    %s\n"%args.output)
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
    if args.output == 'stdout':
        from sys import stdout as outfile
    elif args.output.lower().endswith('.gz'):
        from gzip import open as gopen; outfile = gopen(args.output,'w')
    else:
        outfile = open(args.output,'w')

    # run workflow
    if VERBOSE:
        stderr.write("=== RUN WORKFLOW ===\n"); START = time()
    infile_lines = infile.read().splitlines()
    if infile_lines[0][0] != '>': # Sequences (FASTA)
        raise ValueError("Input file not FASTA")
    if VERBOSE:
        stderr.write("Input file parsed as sequences (FASTA)\n")
    from tempfile import NamedTemporaryFile
    tmp = NamedTemporaryFile(mode='w')
    for l in infile_lines:
        tmp.write(l.strip()); tmp.write('\n')
    tmp.flush()
    if VERBOSE:
        stderr.write("Running MAFFT... "); stderr.flush(); START_MAFFT = time()
    o = run([args.mafft_path, '--auto', '--thread', str(args.num_threads), tmp.name], stdout=PIPE, stderr=PIPE)
    f = open('%s.mafft.aln'%args.input,'wb'); f.write(o.stdout); f.close()
    f = open('%s.mafft.log'%args.input,'wb'); f.write(o.stderr); f.close()
    if VERBOSE:
        stderr.write("%f seconds\n" % (time()-START_MAFFT))
        stderr.write("Running IQ-TREE... "); stderr.flush(); START_IQ = time()
    run([args.iqtree_path, '-s', '%s.mafft.aln'%args.input, '-pre', '%s.iqtree'%args.input, '-nt', str(args.num_threads), '-ntmax', str(args.num_threads), '-redo'], stdout=PIPE, stderr=PIPE)
    if VERBOSE:
        stderr.write("%f seconds\n" % (time()-START_IQ))
        stderr.write("Computing pairwise distances...\n"); stderr.flush()
    dm = read_tree_newick('%s.iqtree.treefile'%args.input).distance_matrix()
    nodes = list(dm.keys())
    for i in range(len(nodes)-1):
        for j in range(i+1,len(nodes)):
            outfile.write('%s,%s,%s\n' % (nodes[i],nodes[j],str(dm[nodes[i]][nodes[j]])))
    outfile.close()
    if VERBOSE:
        END = time(); stderr.write("Total runtime: %s seconds\n" % (END-START))
