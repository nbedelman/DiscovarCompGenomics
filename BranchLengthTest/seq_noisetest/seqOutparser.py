import numpy as numpy
import sys as sys
import ast as ast
import glob
import os
import matplotlib.pyplot as plt

path=sys.argv[1]
estimateArr=[[],[],[]]
for filen in glob.glob(os.path.join(path,'*')):
	#print filen
	F =  open(filen,'r')
	f1=F.readline()
	f1=ast.literal_eval(f1)
	for ii,x in enumerate(f1):
		xi=list(x[0][0])
		xi.sort()
		f1[ii][0][0]=xi
	f1.sort(key=lambda x: x[0][0][0])
	for ii,x in enumerate(f1):
		if x[0][0]==['1','2','3']:
			for z in range(0,3):
				estimateArr[int(x[z][1])-1].append(x[z][3][1])
print estimateArr
plt.hist(estimateArr[1],bins=15)
plt.xlim(0.5,1)
#plt.xlim(2,4)
plt.vlines(0.783,0,15)
#plt.vlines(3,0,17)
#plt.xlabel('Inferred C Value (in coalescent units)')
plt.ylabel('Frequency')
plt.xlabel('Inferred Introgression Proportion')
plt.show()
