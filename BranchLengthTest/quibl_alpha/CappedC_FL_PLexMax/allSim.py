import numpy as numpy
from scipy.stats import expon,truncexpon
import pickle as pickle
import sys as sys

def BLDGenerator(p,samplesize,lmbd,C):
	#Going to do a fixed species tree with given probabilities of introgression.
	outputList=[]
	for x in range (0,samplesize):
		if numpy.random.random_sample()<p:
			outputList.append(expon.rvs()*lmbd)
		else:
			outputList.append(expon.rvs()*lmbd+(C-truncexpon.rvs(C))*lmbd)
#			print (expon.rvs()*lmbd+(C-truncexpon.rvs(C))*lmbd)
	return [[('testA','testB','testC'),outputList,[],[]]]

def BLD_3_Generator(p,q,samplesize,lmbd,C,C2):
	#Going to do a fixed species tree with given probabilities of introgression.
	outputList=[]
	for x in range (0,samplesize):
		if numpy.random.random_sample()<p:
			outputList.append(expon.rvs()*lmbd)
		elif numpy.random.random_sample()<(p+q):
			outputList.append(expon.rvs()*lmbd+(C-truncexpon.rvs(C))*lmbd)
		else:
			outputList.append(expon.rvs()*lmbd+(C2-truncexpon.rvs(C2))*lmbd)
			#print (expon.rvs()*lmbd+(C-truncexpon.rvs(C))*lmbd)
	return [[('testA','testB','testC'),outputList,[],[]]]
#totalOut=BLDGenerator(float(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
#print truncexpon.rvs(10)
#print totalOut

totalOut=BLD_3_Generator(float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]))
pickle.dump(totalOut,open("3disttest.p",'wb'))
