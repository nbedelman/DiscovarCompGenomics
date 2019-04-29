#Gene Tree Introgression Tester
#by Michael Miyagi

import numpy as np
import networkx as nx
import sys
from ete3 import Tree
import glob
import re
import os
import itertools as itt
import cPickle as pickle

def cleanComments(filepath):
#cleanComments removes comments from the newick files at filepath.
	for item in glob.glob(filepath):
		with open(item, 'r+') as tf:
			text=tf.read()
			temp=re.sub("[\[].*?[\]]","",text)
			tf.seek(0)
			tf.write(temp)
			tf.truncate()
			tf.close()

def multi_new(filepath):
#multi_new allows the script to take in files with multiple newick trees within the same file.
	returnlist=[]
	f= open(filepath).read()
	composite=f.split(';')
	for item in composite:
		if item.strip():
			returnlist.append((item+';').strip('\n'))
	return returnlist
def genTripSet(tree):
	leaves=tree.get_leaf_names()
	return itt.combinations(leaves,3)

def inputSet(filepath):
#inputSet takes in a filepath and grabs trees from the input and runs the algorithm, returning a graph.
	print glob.glob(filepath)
	treelist=[]
	t=Tree(multi_new(glob.glob(filepath)[0])[0])
	for fileitem in glob.glob(filepath):
		for item in multi_new(fileitem):
			#t=Tree(item)
			treelist.append(Tree(item))
	return treelist

def iterTrees(tList):
#	print tree
#	leafset=tList[0].get_leaf_names()
	triperate=genTripSet(tList[0])
	lenIt=sum(1 for _ in triperate)
	triperate=genTripSet(tList[0])
	triples=[]
	dist=0
	print lenIt
	for triplet in triperate:
		triples.append(triplet)
	output=np.zeros((4,lenIt))
	output=[[None,[],[],[]] for i in range(lenIt)]
	for tree in tList:
		for index,triplet in enumerate(triples):
	#	print index
	#	print triplet
	#		dist=0 ##HACKY!!
			output[index][0]=triplet
			tempTree=tree.copy('newick')
			tempTree.prune(triplet,preserve_branch_length=True)
			common=tempTree.get_leaves()[0].get_common_ancestor(tempTree.get_leaves())
			for x in common.get_children():
				if not (x.is_leaf()):
					dist=x.get_distance(common)
				else:
					tempInd=triplet.index(x.get_leaf_names()[0])
			output[index][tempInd+1].append(dist)
		
	return output

#t=Tree('((((H,K),(F,I)G),E),((L,(N,Q)O),(P,S)));',format=1)
#print t
#leaves= t.get_leaf_names()
#print leaves
#itt.combinations(leaves,3)
#temp=t.copy('newick')
#temp.prune(('H','K','L'),preserve_branch_length=True)
#print t
#print temp
#print temp.get_topology_id()
#common=temp.get_leaves()[0].get_common_ancestor(temp.get_leaves())
#for x in common.get_children():
#	if not (x.is_leaf()):
#		print x.get_distance(common)

#print #t.get_leaf_names()
#treeList=[]
#treeList.append(t)
f=open(str(sys.argv[2])+'.txt','w')
output=iterTrees(inputSet(sys.argv[1]))
f.write(str(output))
f.close()
pickle.dump(output,open(str(sys.argv[2])+".p","wb"))
