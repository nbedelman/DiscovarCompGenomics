from numpy import *
import sys as sys
import cPickle as pickle
import copy
from joblib import Parallel, delayed
import multiprocessing


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

def qCalc(x,cArr,pArr,lmbd):
	y=[x]
	temp=[]
	for c in range(0,len(cArr)):
		#There was a pArr[c]*pdf...
		temp.append(pArr[c]*pdf(x,cArr[c],lmbd))
	temp=temp/sum(temp)
	y.extend(temp[:])
	return y

def LoglikCalc(XQ,cArr,pArr,lmbd):
	L=0
	for i in range(0,len(pArr)):
		L+=pArr[i]*pdf(XQ[0],cArr[i],lmbd)
	#	print tup[0],cArr[CIndex]
		#if pdf(tup[0],cArr[CIndex],lmbd)<0:
		#	print 'YEESH'
		#	print pdf(tup[0],cArr[CIndex],lmbd),tup[0],cArr[CIndex],lmbd
		#L+=log(prop*(pdf(tup[0],cArr[CIndex],lmbd)))
	return log(L)
def TemplikCalc(XQList,Xdat,XIndex,qIndex,lmbd):
	L=0
	for tup in XQList:
		#MOVED LOG
		L+=log(tup[qIndex+1]*(pdf(tup[0],Xdat[XIndex],lmbd)))
	return L

def cFinder(Xdat,XQList,cArr,lmbd):
	nuC=[0]
	for i in range(1,len(cArr)):
		nuC.append(cSearch(Xdat,XQList,i,lmbd))
	return nuC

def cSearch(Xdat,XQList,ind,lmbd):
	tempLik=[]
	for j in range(0,len(Xdat)):
		tempLik.append(TemplikCalc(XQList,Xdat,j,ind,lmbd))
	op=Xdat[argmax(tempLik)]/lmbd
	return op

def XQListCalc(Xdat,cArr,pArr,lmbd):
	XQList=[]
	for x in Xdat:
		XQList.append(qCalc(x,cArr,pArr,lmbd))
	return XQList

def fullLikCalc(XQList,Xdat,cArr,pArr,lmbd):
#CHANGED BECAUSE OF LIK
	tempSum=0
	for XQ in XQList:
		tempSum+=(LoglikCalc(XQ,cArr,pArr,lmbd))
##		tempsum=(max(tempsum,currSum)+log(exp(tempsum-max(tempsum,currSum))+exp(currSum-max(tempsum,currSum))))
#		tempsum=logaddexp(tempsum,currSum)
	return tempSum

def pUpdate(XQList):
	newpArr=[]
	for i in range(1,len(XQList[0])):
		temp=[]
		for x in XQList:
			temp.append(x[i])
		newpArr.append(mean(temp))
	return newpArr

def monoLUpdate(Xdat,XQList,cArr,pArr,lmbd):
	distro=argmax(pArr)
	submass=0.
	supmass=0.
	for ind,val in enumerate(Xdat):
		if val<cArr[distro]*lmbd:
	#		print val
			submass+=XQList[ind][distro+1]
		else:
			supmass+=XQList[ind][distro+1]
	supmass=len(XQList)-submass
	if submass>0:
		return cArr[distro]*lmbd/(-1*log(1-2*submass/(submass+supmass)))#,[0,(-1*log(1-2*submass/(submass+supmass)))]
	else:
		weighting=[x[0] for x in XQList]
		return average(Xdat,weights=weighting)#,cArr/lmbd*average(Xdat,weights=weighting)

def calcDeriv(XQList,cArr,pArr,lmbd):
##	deriv=0
##	#print XQList
##	for Cindex,C in enumerate(cArr):
##		for x in XQList:
##			if x[0]<C:
##				deriv+=float(x[Cindex+1])*(-1*(x[0]+x[0]*exp(2*x[0]/lmbd)-lmbd+lmbd*exp(2*x[0]/lmbd))/((-1+exp(2*x[0]/lmbd))*pow(lmbd,2)))
##			if x[0]>=C:
##				deriv+=float(x[Cindex+1])*((x[0]-lmbd)/pow(lmbd,2))
##	return deriv
#Redone to fix magnitude errors, assuming 2 distributions or fewer.
	deriv=0
	if len(cArr)==2:
		alpha=cArr[0]*lmbd
		beta=cArr[1]*lmbd
		for x in XQList:
			if x[0]<beta:
				p=pArr[1]
				q=pArr[0]
				deriv+=(x[0]+(beta/(-1+exp(beta/lmbd)))-(exp(2*x[0]/lmbd)*p*(2*x[0]-beta)+(p+2*q)*beta)/((-1+exp(2*x[0]/lmbd))*p+2*(-1+exp(beta/lmbd))*q)-lmbd)/pow(lmbd,2)
			if x[0]>=beta:
				p=pArr[1]
				q=pArr[0]
				deriv+=(x[0]-beta+((p+2*q)*beta)/(p+exp(beta/lmbd)*p+2*q)-lmbd)/pow(lmbd,2)
		return deriv
	elif len(cArr)==1:
		alpha=cArr[0]*lmbd
		for x in XQList:
#			if x[0]<alpha:
#				deriv+=(alpha+(alpha/(-1+exp(alpha/lmbd)))-lmbd-x[0]*coth(x[0]/lmbd))/pow(lmbd,2)
#			if x[0]>=alpha:
#				deriv+=(x+(-1+1/(1+exp(alpha/lmbd)))*alpha-lmbd)/pow(lmbd,2)
			deriv+=(x[0]-lmbd)/pow(lmbd,2)
		return deriv			
	else:
		return null	
	


def badSearch(Xdat,XQList,cArr,pArr,lmbd):
#	submass=0.
#	supmass=0.
#	for ind,val in enumerate(Xdat):
#		if val<cArr[1]*lmbd:
#			submass+=XQList[ind][distro]
#		else:
#			supmass+=XQList[ind][distro]
#	supmass=len(XQList)-submass
#	newlmbd=cArr[distro]*lmbd/(-1*log(1-2*submass/(submass+supmass))) 
	newlmbd=lmbd
	cArr=[x*(lmbd/newlmbd) for x in cArr]
#	print cArr
	likArray=[0]
	likArray.append(fullLikCalc(XQList,Xdat,cArr,pArr,newlmbd))
	k=1
	while (likArray[k]-likArray[k-1])>=0.1*likArray[k] and k<10:
		#print 'im working!'
		biglmbd=newlmbd+newlmbd/2
		smalllmbd=newlmbd/2
		bigcArr=[x*(newlmbd/biglmbd) for x in cArr]
		smallcArr=[x*(newlmbd/smalllmbd) for x in cArr]
		big=fullLikCalc(XQList,Xdat,bigcArr,pArr,biglmbd)
		small=fullLikCalc(XQList,Xdat,smallcArr,pArr,smalllmbd)
		if big > likArray[k] or small > likArray[k]:
			if big>small:
				newlmbd=biglmbd
				cArr=bigcArr
				likArray.append(big)
			elif small>big:
				newlmbd=smalllmbd
				cArr=smallcArr
				likArray.append(small)
		else:
			likArray.append(fullLikCalc(XQList,Xdat,cArr,pArr,newlmbd))
		k+=1
#	print newlmbd
	return newlmbd	
def gradAscent(Xdat,XQList,cArr,pArr,lmbd):
	likArray=[0]
	likArray.append(fullLikCalc(XQList,Xdat,cArr,pArr,lmbd))
	k=1
	stepScale=0.2
	h=1
	while (likArray[k]-likArray[k-1])>=0.05*likArray[k] and k<50 and h<50:
		#print 'im working!'
		newlmbd=lmbd+calcDeriv(XQList,cArr,pArr,lmbd)*stepScale
		if newlmbd<=0:
			stepScale=stepScale/8
			h+=3
		else:
			newcArr=[x*(lmbd/newlmbd) for x in cArr]
			newLik=fullLikCalc(XQList,Xdat,newcArr,pArr,newlmbd)
			if newLik > likArray[k]:
				lmbd=newlmbd
				cArr=newcArr
				likArray.append(newLik)
				k+=1
			else:
				stepScale=stepScale/2
				h+=1
	return lmbd
		

def InitEM(datafp,numDist,thresh,numFlag,data_flag=None):
	#print numDist
	if data_flag is None:
		data=parser(datafp,numFlag)
		data=filter(lambda a: a!=0, data)
	else:
		data=datafp
	cArr=[0]
	pArr=[1./numDist]
##	lmbd=mean(data)
	lmbd=0.02
	for x in range(1,numDist):
#		cArr.append(percentile(data,(x)/float(numDist)*100.)*lmbd)
		cArr.append(lmbd)
		pArr.append(1./numDist)
	XQList=XQListCalc(data,cArr,pArr,lmbd)
	likArray=[0,fullLikCalc(XQList,data,cArr,pArr,lmbd)]
	k=1
	steps=0
	FULLcArr=[cArr]
	FULLpArr=[pArr]
	FULLlmbd=[lmbd]
	while abs(likArray[k]-likArray[k-1])>=thresh*likArray[k] and steps<20:
		steps+=1
		XQList=XQListCalc(data,cArr,pArr,lmbd)
		pArr=pUpdate(XQList)
		FULLpArr.append(pArr)
		cArr=cFinder(data,XQList,cArr,lmbd)
		FULLcArr.append(cArr)
		FULLlmbd.append(lmbd)
		k=k+1
		#likArray.append(fullLikCalc(XQList,data,cArr,pArr,lmbd))
		#pArr=pUpdate(XQList)
		oldlmbd=lmbd
#		lmbd=monoLUpdate(data,XQList,cArr,pArr,lmbd)
#		lmbd=badSearch(data,XQList,cArr,pArr,lmbd)
		lmbd=gradAscent(data,XQList,cArr,pArr,lmbd)
		cArr=[x*(oldlmbd/lmbd) for x in cArr]
		#pArr=pUpdate(XQList)
		likArray.append(fullLikCalc(XQList,data,cArr,pArr,lmbd))	
	index=argmax(likArray)
	return FULLcArr[index-1],FULLpArr[index-1],likArray[index],BIC(likArray[index],len(data),(2*numDist-1)),FULLlmbd[index-1]

#	print cArr,pArr,likArray[-1], BIC(likArray[-1],len(data),(2*(numDist-1)+1)),lmbd
#	return cArr,pArr,likArray[-1], BIC(likArray[-1],len(data),(2*(numDist-1)+1)),lmbd

def parser(filepath,num):
	tempdata=pickle.load(open(filepath,"rb"))
	data=tempdata[num][1]
	return data

def realDataAuto(datapath,outpath,numDist,scaling):
#Given the format of the data: [('a','b','c'),[...],[...],[...]],...
#We want to test every topology from every triplet and output the results.
	print numDist
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
					scaled=[(y)*scaling for y in vect[x] if y>0]
					#print mean(scaled)
					#swap out scaled for vect[x]
					cArr,pArr,MLEst,BIC=InitEM(scaled,numDist,0.001,0,"FULL")
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

def PLDataAuto(datapath,outpath,numDist):
	data=pickle.load(open(datapath,"rb"))
	f=open(outpath,'w')
	num_cores=multiprocessing.cpu_count()
	output=Parallel(n_jobs=num_cores)(delayed(toPar)(ind,vect,numDist) for ind,vect in enumerate(data))
	f.write(str(output))
			

def toPar(ind,vect,numDist):
	out=[]
	fullout=[]
	if vect and ind>=0:
		for x in range(1,4):
			out=[]
			trip=vect[0]
			outgr=vect[0][x-1]
			if vect[x]:
				#scaled=[(y)*float(scaling) for y in vect[x] if y>0]
				cArr,pArr,MLEst,BIC,lmbdEst=InitEM(vect[x],numDist,0.1,0,"FULL")
			else:
				cArr,pArr,MLEst,BIC,lmbdEst=[],[],[],[],[]
			out=[trip,outgr,cArr,pArr,MLEst,BIC,lmbdEst]
			fullout.append(out)
			cArr=[]
			pArr=[]
			MLEst=[]
			BIC=[]
			lmbdEst=[]

	return fullout

def BIC(likelihood, numData, numParams):
	return log(numData)*numParams-2*likelihood

#print InitEM(sys.argv[1],int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))
#accTest2(sys.argv[1],int(sys.argv[2]))
#input file, output file, numdists
PLDataAuto(sys.argv[1],sys.argv[2],int(sys.argv[3]))

