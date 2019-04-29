import numpy as numpy
import sys as sys
import ast as ast

#Format is:
#[triplet #][leaf #][0=triplet, 1=leaf, 2=[cArr], 3=[pArr], 4=logLik, 5=BIC value]

F =  open(sys.argv[1]+"_1Dist.txt",'r')
f1=F.readline()
f1=ast.literal_eval(f1)

J =  open(sys.argv[1]+"_2Dist.txt",'r')
f2=J.readline()
f2=ast.literal_eval(f2)

for ii,x in enumerate(f1):
	x[0][0]=list(x[0][0])
	x[0][0].sort()
	f2[ii][0][0]=list(f2[ii][0][0])
	f2[ii][0][0].sort()
	x.sort(key=lambda y: y[1])
	f2[ii].sort(key=lambda y: y[1])

f1.sort(key=lambda x: x[0][0])
f2.sort(key=lambda x: x[0][0])

#To Check C value concordance.
f1Par=[]
f2Par=[]
for ii,x in enumerate(f1):
	if x[0][2] and f2[ii][0][2] :
		f1Par.append(x[0][2][0])
		f2Par.append(f2[ii][0][2][0])
	if x[1][2] and f2[ii][1][2] :
		f1Par.append(x[1][2][0])
		f2Par.append(f2[ii][1][2][0])
	if x[2][2] and f2[ii][2][2] :
		f1Par.append(x[2][2][0])
		f2Par.append(f2[ii][2][2][0])
#print f1Par
#print f2Par
#print numpy.corrcoef(f1Par,f2Par)

#To check BIC assignment.
for ii,x in enumerate(f1):
	print x[0][0]
	print x[0][2],x[1][2],x[2][2]
	print f2[ii][0][2],f2[ii][1][2],f2[ii][2][2]
	if x[0][5] > f2[ii][0][5]:
		print True
	else: print False
	if x[1][5] > f2[ii][1][5]:
		print True
	else: print False
	if x[2][5] > f2[ii][2][5]:
		print True
	else: print False

def treePainter():
#In principle- given a full tree and an input distribution, paint the edges by distribution


def chomPainter():
#Given a series of triplet trees, Paint a linear chomosome by distribution of that the trees belong to.




#Check if the triplet order is preserved.
#for ii,x in enumerate(f1):
#	print list(x[0][0])
#	print f2[ii][0][0]
#	print (x[0][0])==(f2[ii][0][0])
#Answer: Yes, after sorting!
