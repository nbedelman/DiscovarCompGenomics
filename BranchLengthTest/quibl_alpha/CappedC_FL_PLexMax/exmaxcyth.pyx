import numpy as np
from libc.math cimport exp,sinh,cosh,log,pow
cimport numpy as np
DTYPE=np.double

def pdf(double x,double C,double lmbd):
	cdef double y
	if x>0 and x>=C*lmbd:
		y=0.5/lmbd*exp(-x/lmbd)*(1+exp(C))
	elif x>0 and C*lmbd>x:
		y=sinh(x/lmbd)/(lmbd*(-1+exp(C)))
	else:
		y=0
	return y
def pUpdate(list XQList):
	newpArr=[]
	for i in range(1,len(XQList[0])):
		temp=[]
		for x in XQList:
			temp.append(x[i])
		newpArr.append(np.mean(temp))
	return newpArr

def TemplikCalc(list XQList, double testC, int CIndex, double lmbd, list cArr):#, dict lookup):
	cdef double L=0
	cdef double temp=0
	cdef double y=0
	cdef double x=0
	cdef double C=0
	cdef int ind,indi
	cdef double pival=0
	for ind,tup in enumerate(XQList):
		temp=0
		#for indi,pival in enumerate(tup[1:]):
		for indi in range(len(tup[1:])):
			#print indi
			pival=tup[1:][indi]#[indi+1]
			if indi==CIndex:
				x=tup[0]
				C=testC
				if x>0 and x>=C*lmbd:
					y=0.5/lmbd*exp(-x/lmbd)*(1+exp(C))
				elif x>0 and C*lmbd>x:
					y=sinh(x/lmbd)/(lmbd*(-1+exp(C)))
				else:
					y=0
				temp+=pival*y
				#temp=temp+pival*(pdf(tup[0],testC,lmbd))
			else:
				x=tup[0]
				C=cArr[indi]
				if x>0 and x>=C*lmbd:
					y=0.5/lmbd*exp(-x/lmbd)*(1+exp(C))
				elif x>0 and C*lmbd>x:
					y=sinh(x/lmbd)/(lmbd*(-1+exp(C)))
				else:
					y=0
				temp+=pival*y
				#temp+=pival*(pdf(tup[0],cArr[indi],lmbd))
				#temp+=lookup[str(ind)+str(indi)]
		L+=log(temp)
	return L

def cLookup(list XQList, list cArr, double lmbd):
	cdef dict lookup={}
	for ind,tup in enumerate(XQList):
		for indi,pival in enumerate(tup[1:]):
			lookup[str(ind)+str(indi)]=pival*(pdf(tup[0],cArr[indi],lmbd))
	return lookup

def cSearch(double[:] Xdat,list XQList,int ind,double lmbd,list cArr, float threshold):
	#RESCALEXdat=[x/lmbd for x in Xdat]
	cdef int datlen=Xdat.size
	cdef np.ndarray RESCALEXdat=np.zeros(datlen,dtype=float)
	cdef double op=0
	cdef int x,j
	for x in range(datlen):
		RESCALEXdat[x]=Xdat[x]/lmbd
	#tempLik=[]
	cdef np.ndarray tempLik=np.zeros(datlen)
	#lookup=cLookup(XQList,cArr,lmbd)
	for j in range(0,datlen):
		if RESCALEXdat[j]>threshold:
			tempLik[j]=TemplikCalc(XQList,RESCALEXdat[j],ind,lmbd,cArr)#,lookup))
		else:
			tempLik[j]==-np.NINF
	op=RESCALEXdat[np.argmax(tempLik)]
	return op

def cFinder(double[:] Xdat,list XQList,list cArr,double lmbd, float threshold):
	cdef int clen=len(cArr)
	cdef list nuC=[0]
	#print clen
	for i in range(1,clen):
		nuC.append(cSearch(Xdat,XQList,i,lmbd,cArr,threshold))
	return nuC
def qCalc(double x, list cArr,list pArr,double lmbd):
	y=[x]
	temp=[]
	for c in range(0,len(cArr)):
		temp.append(pArr[c]*pdf(x,cArr[c],lmbd))
	summand=sum(temp)
	temp=[z/summand for z in temp]
	y.extend(temp[:])
	return y
def XQListCalc(double[:] Xdat,list cArr,list pArr,double lmbd):
	XQList=[]
	for x in Xdat:
		XQList.append(qCalc(x,cArr,pArr,lmbd))
	return XQList
def LoglikCalc(double X,list cArr,list pArr,double lmbd):
	cdef double L
	L=0
	for i in range(0,len(pArr)):
		#print pArr,cArr
		L+=pArr[i]*pdf(X,cArr[i],lmbd)
	return log(L)

def fullLikCalc(list XQList,double[:] Xdat,list cArr,list pArr,double lmbd):
	tempSum=0
	for X in Xdat:
		tempSum+=(LoglikCalc(X,cArr,pArr,lmbd))
	return tempSum

def calcDeriv(list XQList,list cArr,list pArr,double lmbd):
	cdef double val
	cdef double alpha
	cdef double beta
	deriv=0
	if len(cArr)==2:
		alpha=0#cArr[0]*lmbd
		beta=cArr[1]*lmbd
		for x in XQList:
			val=x[0]
			if val<beta:
				p=x[2]#pArr[1] #CHECK THIS?
				q=x[1]#pArr[0]
				deriv+=(val+(beta/(-1+exp(beta/lmbd)))-(exp(2*val/lmbd)*p*(2*val-beta)+(p+2*q)*beta)/((-1+exp(2*val/lmbd))*p+2*(-1+exp(beta/lmbd))*q)-lmbd)/pow(lmbd,2)
			if val>=beta:
				p=x[2]#pArr[1]
				q=x[1]#pArr[0]
				deriv+=(x[0]-beta+((p+2*q)*beta)/(p+exp(beta/lmbd)*p+2*q)-lmbd)/pow(lmbd,2)
		return deriv
	elif len(cArr)==1:
		alpha=cArr[0]*lmbd
		for x in XQList:
			deriv+=(x[0]-lmbd)/pow(lmbd,2)
		return deriv			
	elif len(cArr)==3:
		b=cArr[1]*lmbd
		g=cArr[2]*lmbd
		for x in XQList:
			p=x[2]#pArr[1]
			q=x[3]#pArr[0]
			if x[0]<b and x[0]<g:
				deriv+=(lmbd*(-2+ p + q) - (b + lmbd)*p*exp(b/lmbd) + x[0]*(2 + p*(-1 + exp(b/lmbd)) + q*(-1 + exp(g/lmbd))) - (g + lmbd)*q*exp(g/lmbd))/(pow(lmbd,2)*(2 + p*(-1 + exp(b/lmbd)) + q*(-1 + exp(g/lmbd))))
			elif x[0]>=b and x[0]<g:
				deriv+=(-((q*x[0]*cosh(x[0]/lmbd))/(lmbd**3*(-1 + exp(g/lmbd)))) - ((1 - p - q)*exp(-(x[0]/lmbd)))/lmbd**2 + ((1 - p - q)*x[0]*exp(-(x[0]/lmbd)))/lmbd**3 - (p*(1 + exp(b/lmbd))*exp(-(x[0]/lmbd)))/(2.*lmbd**2) + (p*x[0]*(1 + exp(b/lmbd))*exp(-(x[0]/lmbd)))/(2.*lmbd**3) - (b*p*exp(b/lmbd - x[0]/lmbd))/(2.*lmbd**3) - (q*sinh(x[0]/lmbd))/(lmbd**2*(-1 + exp(g/lmbd))) + (g*q*exp(g/lmbd)*sinh(x[0]/lmbd))/(lmbd**3*(-1 + exp(g/lmbd))**2))/(((1 - p - q)*exp(-(x[0]/lmbd)))/lmbd + (p*(1 + exp(b/lmbd))*exp(-(x[0]/lmbd)))/(2.*lmbd) + (q*sinh(x[0]/lmbd))/(lmbd*(-1 + exp(g/lmbd))))
			elif x[0]>=b and x[0]>=g:
				deriv+=(x[0]*cosh(x[0]/lmbd)*(-(p/(-1 + exp(b/lmbd))) - q/(-1 + exp(g/lmbd))) - (-1 + p + q)*(-lmbd + x[0])*exp(-(x[0]/lmbd)) + ((g*q*1/sinh(g/(2.*lmbd))**2)/4. + (p*(lmbd + (b - lmbd)*exp(b/lmbd)))/(-1 + exp(b/lmbd))**2 - (lmbd*q)/(-1 + exp(g/lmbd)))*sinh(x[0]/lmbd))/(lmbd**2*(-((-1 + p + q)*exp(-(x[0]/lmbd))) + (p/(-1 + exp(b/lmbd)) + q/(-1 + exp(g/lmbd)))*sinh(x[0]/lmbd)))
			elif x[0]>=g and x[0]<b:
				temp=g
				g=b
				b=temp
				p=x[3]
				q=x[2]
				deriv+=(-((q*x[0]*cosh(x[0]/lmbd))/(lmbd**3*(-1 + exp(g/lmbd)))) - ((1 - p - q)*exp(-(x[0]/lmbd)))/lmbd**2 + ((1 - p - q)*x[0]*exp(-(x[0]/lmbd)))/lmbd**3 - (p*(1 + exp(b/lmbd))*exp(-(x[0]/lmbd)))/(2.*lmbd**2) + (p*x[0]*(1 + exp(b/lmbd))*exp(-(x[0]/lmbd)))/(2.*lmbd**3) - (b*p*exp(b/lmbd - x[0]/lmbd))/(2.*lmbd**3) - (q*sinh(x[0]/lmbd))/(lmbd**2*(-1 + exp(g/lmbd))) + (g*q*exp(g/lmbd)*sinh(x[0]/lmbd))/(lmbd**3*(-1 + exp(g/lmbd))**2))/(((1 - p - q)*exp(-(x[0]/lmbd)))/lmbd + (p*(1 + exp(b/lmbd))*exp(-(x[0]/lmbd)))/(2.*lmbd) + (q*sinh(x[0]/lmbd))/(lmbd*(-1 + exp(g/lmbd))))	
		return deriv
	else:
		return None	

def gradAscent(np.ndarray Xdat,list XQList,list cArr,list pArr,double lmbd, double SCALAR):
	if len(cArr)==1:
		return np.mean(Xdat),SCALAR
	elif len(cArr)>1:
		likArray=[0]
		likArray.append(fullLikCalc(XQList,Xdat,cArr,pArr,lmbd))
		k=1
		stepScale=SCALAR
		h=1
		while abs(likArray[k]-likArray[k-1])>=0.001*stepScale*abs(likArray[k]) and k<200 and h<200:
			newlmbd=lmbd+calcDeriv(XQList,cArr,pArr,lmbd)*stepScale
			if newlmbd<=0:
				stepScale=stepScale/2
				h+=1
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
		return lmbd,stepScale

