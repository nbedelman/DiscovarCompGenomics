#!/usr/bin/env python
import sys
from StringIO import StringIO
from Bio import Phylo
from Bio.Phylo.Consensus import _BitString

#for full erato clade

treeFile=sys.argv[1]
outFile=sys.argv[2]


#the input here is a single columns of newick trees.
#The output will be a csv with two columns: region and type.

# function to convert a list of nwk trees into a ETE tree objects
def makeTree(tree):
    tree = tree.rstrip()
    if not tree == 'NA':
    #     tree = Tree(tree)
        tree=Phylo.read(StringIO(tree), 'newick')
        return (tree)
    else:
        return ('NA')

def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs

def compare(tree1, tree2):
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False


t=open(treeFile, "r")
out=open(outFile, "a+")

foundTrees=[]
foundTreeNum=0
foundTreeCount={}

for entry in t:
    tree=makeTree(entry.strip())
    for i in tree.get_terminals():
        if i.name=='Etal':
            Etal=i
    tree.root_with_outgroup(Etal)
    matchName=''
    for i in foundTrees:
        if compare(tree, i):
            matchName=i.name
            foundTreeCount[i.name][1]+=1
            out.write('''%s,%s\n''' % (region,matchName))
            break
        if matchName=='':
            tree.name="tree"+str(foundTreeNum)
            foundTrees.append(tree)
            foundTreeCount[tree.name]=[tree,1]
            foundTreeNum+=1
            out.write('''%s,%s\n''' % (region,tree.name))
    else:
    	out.write('''%s,%s\n''' % (region,"NA"))

f=open(outFile+".trees", "a+")

for i in foundTreeCount.keys():
    f.write('''%s\t%i\t%s''' % (i, foundTreeCount[i][1],foundTreeCount[i][0].format("newick")))

f.close()
out.close()
t.close()
