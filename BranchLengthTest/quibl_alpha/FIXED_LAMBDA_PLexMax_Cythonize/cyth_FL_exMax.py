#Michael Miyagi
#8/21 
#Cythonizing

from numpy import *
import sys as sys
import cPickle as pickle
import copy
from joblib import Parallel, delayed
import multiprocessing
import cProfile as cProfile
import exmaxcyth as emc
import pstats,StringIO

def pdf(x,C,lmbd):
	if x>0 and x>=C*lmbd:
		y=0.5/lmbd*exp(-x/lmbd)*(1+exp(C))
		#print y,lmbd,x,C
	elif x>0 and C*lmbd>x:
		y=sinh(x/lmbd)/(lmbd*(-1+exp(C)))
		#print y
	else:
		y=0
	#print y
	return y

#def qCalc(x,cArr,pArr,lmbd):
#	y=[x]
#	temp=[]
#	for c in range(0,len(cArr)):
#		temp.append(pArr[c]*pdf(x,cArr[c],lmbd))
#	temp=temp/sum(temp)
#	y.extend(temp[:])
#	return y

#def LoglikCalc(X,cArr,pArr,lmbd):
#	L=0
#	for i in range(0,len(pArr)):
#		L+=pArr[i]*pdf(X,cArr[i],lmbd)
#	return log(L)

#def TemplikCalc(XQList,testC,CIndex,lmbd,cArr):
#	L=0
#	for tup in XQList:
#		temp=0
#		for indi,pival in enumerate(tup[1:]):
#			if indi==CIndex:
#				temp+=pival*(pdf(tup[0],testC,lmbd))
#			else:
#				temp+=pival*(pdf(tup[0],cArr[indi],lmbd))
#		L+=log(temp)
#	return L

#def cFinder(Xdat,XQList,cArr,lmbd):
#	nuC=[0]
#	for i in range(1,len(cArr)):
#		nuC.append(emc.cSearch(Xdat,XQList,i,lmbd,cArr))
#	return nuC

#def cSearch(Xdat,XQList,ind,lmbd,cArr):
#	RESCALEXdat=[x/lmbd for x in Xdat]# if x>cArr[ind-1]]
#	tempLik=[]
#	#print RESCALEXdat	
#	for j in range(0,len(RESCALEXdat)):
#		tempLik.append(emc.TemplikCalc(XQList,RESCALEXdat[j],ind,lmbd,cArr))
#	#print tempLik
#	op=RESCALEXdat[argmax(tempLik)]#/lmbd
#	return op

#def XQListCalc(Xdat,cArr,pArr,lmbd):
#	XQList=[]
#	for x in Xdat:
#		XQList.append(qCalc(x,cArr,pArr,lmbd))
#	return XQList

#def fullLikCalc(XQList,Xdat,cArr,pArr,lmbd):
#	tempSum=0
#	for X in Xdat:
#		tempSum+=(LoglikCalc(X,cArr,pArr,lmbd))
#	return tempSum

#def pUpdate(XQList):
#	newpArr=[]
#	for i in range(1,len(XQList[0])):
#		temp=[]
#		for x in XQList:
#			temp.append(x[i])
#		newpArr.append(mean(temp))
#	return newpArr


def InitEM(datafp,numDist,thresh,numSteps,lmbd,numFlag=None,data_flag=None):
	#print numDist
	#if data_flag is None:
	#	data=parser(datafp,numFlag)
	#	data=filter(lambda a: a!=0, data)
	#else:
	data=datafp
	#print type(data)
	cArr=list(zeros(numDist))
	pArr=[1./numDist]
	lmbd=lmbd#mean(data)
	for x in range(1,numDist):
		cArr[x]=x
		pArr.append(1./numDist)
	XQList=emc.XQListCalc(data,cArr,pArr,lmbd)
#	print XQList
#	print cArr
	likArray=[0,emc.fullLikCalc(XQList,data,cArr,pArr,lmbd)]
	k=1
	steps=0
	FULLcArr=[cArr]
	FULLpArr=[pArr]
	FULLlmbd=[lmbd]
	stepScale=0.5
	while steps<numSteps:# and abs(likArray[k]-likArray[k-1])>=thresh*abs(likArray[k]):# and steps<100:
		steps+=1
		XQList=emc.XQListCalc(data,cArr,pArr,lmbd)
		pArr=emc.pUpdate(XQList)
		FULLpArr.append(pArr)
		cArr=emc.cFinder(data,XQList,cArr,lmbd)
#		oldlmbd=lmbd
#		lmbd,stepScale=emc.gradAscent(data,XQList,cArr,pArr,lmbd,stepScale)
#		cArr=[x*(oldlmbd/lmbd) for x in cArr]
		#cArr=multiply(cArr,oldlmbd/lmbd)
		FULLcArr.append(cArr)
		FULLlmbd.append(lmbd)
		likArray.append(emc.fullLikCalc(XQList,data,cArr,pArr,lmbd))
	if likArray:
		#likArray.append(emc.fullLikCalc(XQList,data,[0,2,4],[0.2,0.4,0.4],0.1))
		index=argmax(likArray)
		return FULLcArr[index-1],FULLpArr[index-1],likArray,BIC(likArray[index],len(data),(2*numDist-1)),FULLlmbd[index-1],index
	else:
		return None	

def parser(filepath,num):
	tempdata=pickle.load(open(filepath,"rb"))
	data=tempdata[num][1]
	return data

def realDataAuto(datapath,outpath,numDist,scaling,numSteps):
#Given the format of the data: [('a','b','c'),[...],[...],[...]],...
#We want to test every topology from every triplet and output the results.
#	print numDist
	data=pickle.load(open(datapath,"rb"))
	f=open(outpath,'w')
	fullout=[]
	for ind,vect in enumerate(data):
		if vect and ind>=0:
			for x in range(1,4):
				out=[]
				trip=vect[0]
				outgr=vect[0][x-1]
				if vect[x]:
					#scaled is to test out linear scaling.
					scaled=[y*scaling for y in vect[x] if y>0]
					#print mean(scaled)
					#swap out scaled for vect[x]
					cArr,pArr,MLEst,BIC=InitEM(scaled,numDist,0.05,numSteps)
				else:
					cArr,pArr,MLEst,BIC=[],[],[],[]
				out=[trip,outgr,cArr,pArr,MLEst,BIC]
				f.write(str(out)+'\n')
				f.flush()
				cArr=[]
				pArr=[]
				MLEst=[]
				BIC=[]
				fullout.append(out)
		print ind
#	pickle.dump(fullout,open("FULL_PHYML_exmax.p","wb"))

def PLDataAuto(datapath,outpath,numDist,numSteps,lmbd):
	data=pickle.load(open(datapath,"rb"))
	data=asarray(data)
	f=open(outpath,'w')
	num_cores=multiprocessing.cpu_count()
	output=Parallel(n_jobs=num_cores)(delayed(toPar)(ind,vect,numDist,numSteps,lmbd) for ind,vect in enumerate(data))
	f.write(str(output))

def bootStrapDataAuto(datapath,outpath,numDist):
	data=pickle.load(open(datapath,"rb"))
	f=open(outpath,'w')
	num_cores=multiprocessing.cpu_count()
	output=Parallel(n_jobs=num_cores)(delayed(bootStrapPar)(ind,vect,numDist) for ind,vect in enumerate(data))
	f.write(str(output))

def bootstrapify(data,iterations):
	output=[]
	for x in range(0,iterations):
		output.append(random.choice(data,size=len(data),replace=True))
	return output

def bootStrapPar(ind,vect,numDist,numSteps):
	out=[]
	fullout=[]
	if vect and ind>=0:
		for x in range(1,4):
			out=[]
			trip=vect[0]
			outgr=vect[0][x-1]
			remZero=[y for y in vect[x] if y!=0]
			cVect=[]
			pVect=[]
			mVect=[]
			bVect=[]
			lVect=[]
			iVect=[]
			if len(remZero)>=10: #vect[x]:
				#scaled=[(y)*float(scaling) for y in vect[x] if y>0]
				BSdata=bootstrapify(remZero,10)
				for y in range(0,len(BSdata)):
					temp=InitEM(BSdata[y],numDist,0.1,numSteps)
					cVect.append(temp[0])
					pVect.append(temp[1])
					mVect.append(temp[2])
					bVect.append(temp[3])
					lVect.append(temp[4])
					iVect.append(temp[5])
				#cArr,pArr,MLEst,BIC,lmbdEst,maxIndex=InitEM(remZero,numDist,0.1,0,"FULL")
			else:
				cVect,pVect,mVect,bVect,lVect,iVect=[],[],[],[],[],[]
			out=[trip,outgr,cVect,pVect,mVect,bVect,lVect,len(remZero),iVect]
			fullout.append(out)
			cVect=[]
			pVect=[]
			mVect=[]
			bVect=[]
			lVect=[]
			iVect=[]
	return fullout

def toPar(ind,vect,numDist,numSteps,lmbd):
	out=[]
	fullout=[]
	if vect.size!=0 and ind>=0:
		for x in range(1,4):
			out=[]
			trip=vect[0]
			outgr=vect[0][x-1]
			remZero=[y for y in vect[x] if y!=0]
			remZero=asarray(remZero)
			if remZero.size>=10: #vect[x]:
				#scaled=[(y)*float(scaling) for y in vect[x] if y>0]
				cArr,pArr,MLEst,BIC,lmbdEst,maxIndex=InitEM(remZero,numDist,0.1,numSteps,lmbd)
			else:
				cArr,pArr,MLEst,BIC,lmbdEst,maxIndex=[],[],[],[],[],[]
			out=[trip,outgr,cArr,pArr,MLEst,BIC,lmbdEst,len(remZero),maxIndex]
			fullout.append(out)
			cArr=[]
			pArr=[]
			MLEst=[]
			BIC=[]
			lmbdEst=[]
	return fullout


def MeterDataAuto(datapath,outpath,numDist,numSteps):
	data=pickle.load(open(datapath,"rb"))
	f=open(outpath,'w')
	num_cores=multiprocessing.cpu_count()
	output=[]
	data=asarray(data)
	#print type(data)
	#print type(data[0])
	for ind,vect in enumerate(data):
		output.append(toPar(ind,vect,numDist,numSteps))
	f.write(str(output))



def BIC(likelihood, numData, numParams):
	return log(numData)*numParams-1-2*likelihood

#print InitEM(sys.argv[1],int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))
#accTest2(sys.argv[1],int(sys.argv[2]))
#input file, output file, numdists
#pr=cProfile.Profile()
#pr.enable()
PLDataAuto(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),float(sys.argv[5]))
#pr.disable()
#s=StringIO.StringIO()
#sortby='ncalls'
#ps=pstats.Stats(pr,stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()
#bootStrapDataAuto(sys.argv[1],sys.argv[2],int(sys.argv[3]))

