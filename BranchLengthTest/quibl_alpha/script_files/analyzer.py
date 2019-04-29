from numpy import *
from scipy import stats
import sys as sys
import cPickle as pickle
import matplotlib.pyplot as plt

#data=pickle.load(open(sys.argv[1],"rb"))
#print data[556][0]
#for i in range(1,4):
#	plt.figure(i)
#	plt.hist(data[556][i])
#plt.show()


def maxlenCorr(data,thresh):
	op=[0,0]
	for x in range(1,len(data)):
		if (len(data[x][1])>thresh and len(data[x][2])>thresh and len(data[x][3])>thresh):
			print data[x][0],x
			if max(data[x][1])>max(data[x][2]) and max(data[x][1])>max(data[x][3]) and len(data[x][1])>len(data[x][2]) and len(data[x][1])>len(data[x][3]):
				op[0]+=1
			elif max(data[x][2])>max(data[x][1]) and  max(data[x][2])>max(data[x][3]) and len(data[x][2])>len(data[x][1]) and len(data[x][2])>len(data[x][1]):
				op[0]+=1
			elif max(data[x][3])>max(data[x][2]) and  max(data[x][3])>max(data[x][1]) and len(data[x][3])>len(data[x][2]) and len(data[x][3])>len(data[x][1]):
				op[0]+=1
	
			else: op[1]+=1
	print op

def rankTest(arg):
	ou=[]
	ou.append(stats.kruskal(data[arg][1],data[arg][2],data[arg][3])[1])
	ou.append(stats.mannwhitneyu(data[arg][1],data[arg][2])[1])
	ou.append(stats.mannwhitneyu(data[arg][1],data[arg][3])[1])
	ou.append(stats.mannwhitneyu(data[arg][2],data[arg][3])[1])
	return ou
def fullTester():
	fullou=[]
	counter =0
	for x in range(1,559):
		if (len(data[x][1])>10 and len(data[x][2])>10 and len(data[x][3])>10):
			counter+=1
			if rankTest(x)[0]<0.05:
				fullou.append(rankTest(x))		
	print counter
	return fullou	
def getLengths(data):
	for x in data:
		for y in range(1,4):
			print len(x[y])
	
#x=fullTester()
#print x
#print len(x)
#maxlenCorr(pickle.load(open(sys.argv[1],"rb")),int(sys.argv[2]))
getLengths(pickle.load(open(sys.argv[1],"rb")))
