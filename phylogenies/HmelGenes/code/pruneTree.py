
from ete3 import Tree
import argparse

pruneArgs=argparse.ArgumentParser()
pruneArgs.add_argument("-t","--tree", type=str)
pruneArgs.add_argument("-s","--species",nargs='+')
pruneArgs.add_argument("-o","--outFile", type=str)

args = pruneArgs.parse_args()

from ete3 import Tree
t = Tree(args.tree)
t.prune(args.species, preserve_branch_length=True)
t.write(outfile=args.outFile)
