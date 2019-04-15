#!/usr/bin/env python
import sys
from StringIO import StringIO
from Bio import Phylo
from Bio.Phylo.Consensus import _BitString

#for full erato clade

treeFile=sys.argv[1]
outFile=sys.argv[2]

def makeTree(tree):
    tree = tree.rstrip()
    if not tree == 'NA':
    #     tree = Tree(tree)
        tree=Phylo.read(StringIO(tree), 'newick')
        return (tree)
    else:
        return ('NA')

t=open(treeFile, "r")
out=open(outFile, "a+")

for entry in t:
    atts=entry.split()
    region=entry.split()[0]
    if len(atts)!=3:
        out.write('''%s\t%s\n''' % (region,"NA"))
    else:
        region=entry.split()[0]
        tree=makeTree(atts[1])
        likelihood=atts[2]
        if tree != 'NA':
            for i in tree.get_terminals():
                if i.name=='Etal':
                    Etal=i
            tree.root_with_outgroup(Etal)
            out.write('''%s\t%s\t%s\n''' % (region, tree.format("newick").strip(), likelihood))
t.close()
out.close()
