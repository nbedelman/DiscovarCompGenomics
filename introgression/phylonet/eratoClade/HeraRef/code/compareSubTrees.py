#!/usr/bin/env python

#IT's really hard to figure out how closely related two networks are.
#In phylonet you can use the Charnet command to extract all the bifurcating paths through the network
#Maybe the answer is to see how common each bifurcating path is? Worth a try at least.
#borrowing heavily from findCommonTrees.py and cactusToTrees.py

import sys
from Bio import Phylo

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _Matrix
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from StringIO import StringIO
from Bio.Phylo.Consensus import _BitString
import numpy
import pylab
import sys

charNetFile=sys.argv[1]
outBase=sys.argv[2]

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

def groupTopologies(charNetFile):
    commonTreeDict={}
    networkAttsDict={}
    foundTreeNum=0
    f=open(charNetFile,"r")
    thisNet=''
    for line in f:
        if "nexus" in line:
            thisNet=line.strip()
            networkAttsDict[thisNet]=0
        else:
            tree=makeTree(line)
            for i in tree.get_terminals():
                if i.name=='Etal':
                    Etal=i
            tree.root_with_outgroup(Etal)
            matchName=''
            for commonTree in commonTreeDict.keys():
                if compare(tree, commonTreeDict[commonTree][2]):
                    matchName=commonTree
                    if 'qpGraph' in thisNet:
                        commonTreeDict[commonTree][0]+=1
                    else:
                        commonTreeDict[commonTree][1]+=1
                    break
            if matchName=='':
                if 'qpGraph' in thisNet:
                    commonTreeDict["topology"+str(foundTreeNum)]=[1,0,tree]
                else:
                    commonTreeDict["topology"+str(foundTreeNum)]=[0,1,tree]
                foundTreeNum+=1
            networkAttsDict[thisNet]+=1
    return [commonTreeDict,networkAttsDict]

def writeCommonTreeDict(treeDict,fileName):
    o=open(fileName,"w")
    for key in treeDict:
        o.write('''%s\t%s\t%s\t%s''' % (key, treeDict[key][0],treeDict[key][1],treeDict[key][2].format("newick")))

def writeNetworkAttsDict(treeDict,fileName):
    o=open(fileName,"w")
    for key in treeDict:
        o.write('''%s\t%i\n''' % (key, treeDict[key]))

charDicts=groupTopologies(charNetFile)
commonTreeDict=charDicts[0]
networkAttsDict=charDicts[1]
writeCommonTreeDict(commonTreeDict,outBase+"_commonTrees.tsv")
writeNetworkAttsDict(networkAttsDict, outBase+"_networkAtts.tsv")
