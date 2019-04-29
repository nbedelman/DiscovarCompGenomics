import numpy as np
import networkx as nx
import sys
from ete3 import Tree
import glob
import re
import os
import itertools as itt
import cPickle as pickle
import ast as ast


F =  open(sys.argv[1],'r')
f1=F.readline()
f1=ast.literal_eval(f1)
searchTrip=['HeraRef','Hhsa','Htel']#['1','2','3']

for ii,x1 in enumerate(f1):
	x1=list(x1[0][0])

	x1.sort()

	f1[ii][0][0]=x1
	#print f1[ii][0][0]
lookup=[[],[],[]]
for ii,x in enumerate(f1):
	if x[0][0]==searchTrip:
		lookup[0]=[x[0][1],x[0][2][1],x[0][3][1],x[0][-3]]
		lookup[1]=[x[1][1],x[1][2][1],x[1][3][1],x[0][-3]]
		lookup[2]=[x[2][1],x[2][2][1],x[2][3][1],x[0][-3]]

def pdf(x,C,lmbd):
	if x>0 and x>=C*lmbd:
		y=0.5/lmbd*np.exp(-x/lmbd)*(1+np.exp(C))
		#print y,lmbd,x,C
	elif x>0 and C*lmbd>x:
		y=np.sinh(x/lmbd)/(lmbd*(-1+np.exp(C)))
		#print y
	else:
		y=0
	#print y
	return y


def multi_new(filepath):
#multi_new allows the script to take in files with multiple newick trees within the same file.
	returnlist=[]
	f= open(filepath).read()
	composite=f.split(';')
	for item in composite:
		if item.strip():
			returnlist.append((item+';').strip('\n'))
	return returnlist

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

def genTripSet(tree):
	leaves=tree.get_leaf_names()
	print leaves
	return itt.combinations(leaves,3)


def iterTrees(tList):
	triples=[['HeraRef','Hhsa','Htel']]
	lenIt=1
	#output=np.zeros((4,lenIt))
	#output=[[None,[],[],[]] for i in range(lenIt)]
	output=[0]*len(tList)
	#print len(tList)
	#print output
	for TIND,tree in enumerate(tList):
#		tree.set_outgroup('4')
		for index,triplet in enumerate(triples):
			#output[index][0]=triplet
			tempTree=tree.copy('newick')
			tempTree.prune(triplet,preserve_branch_length=True)
			common=tempTree.get_leaves()[0].get_common_ancestor(tempTree.get_leaves())
			for x in common.get_children():
				if not (x.is_leaf()):
					dist=x.get_distance(common)
				else:
					tempInd=triplet.index(x.get_leaf_names()[0])
			#output[index][tempInd+1].append([TIND,prob(dist,tempInd+1)])
			output[TIND]=[TIND,tempInd+1,dist,prob(tempInd+1,dist)]
	return output

#print lookup
triplE=['HeraRef','Hhsa','Htel']
def prob(tempInd,dist):
	for y in lookup:
		#print y
		#if int(y[0])==int(tempInd):
		if y[0]==triplE[tempInd-1]:
		#	print y[2]*pdf(dist,y[1],y[3])
			#print 'hi'
			return y[2]*pdf(dist,y[1],y[3])/((1-y[2])*pdf(dist,0,y[3])+y[2]*pdf(dist,y[1],y[3]))



#f=open(str(sys.argv[2])+'.txt','w')
output=iterTrees(inputSet(sys.argv[2]))
print output


