from numpy import *
from scipy.stats import expon,iqr
import sys as sys
import cPickle as pickle
import matplotlib.pyplot as plt

data=pickle.load(open(sys.argv[1],"rb"))
param=int(sys.argv[2])
n_bins=200
#print data[param][0]
#print len(data[param][1])
#for i in range(1,4):
#	plt.figure(i)
#	plt.hist(data[param][i])
#plt.show()
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

#param2=int(sys.argv[3])
#stackdat=data[param][param2]#,data[param][2],data[param][3]]
#scalar=1./mean(stackdat)
#print mean(stackdat)*mean(stackdat), var(stackdat)
#seconddat=data[param][3]
#scaledat=[i*scalar for i in seconddat]
#plt.hist(scaledat,normed=0,histtype='bar',stacked=False)
#

for ii,y in enumerate(data):
	x=y[0]
	#print x
	#print index, element
	x=list(x)
	x.sort()
	data[ii][0]=x
#print data[0][0][0]
data.sort(key=lambda x: x[0][0])
#print transpose(data)[0]

for ii,y in enumerate(data):
	print ii,y[0]

print data[param][0]
stackscale1=[x for x in data[param][1]]
stackscale2=[x for x in data[param][2]]
stackscale3=[x for x in data[param][3]]
print len(stackscale1),len(stackscale2),len(stackscale3)

fig,ax=plt.subplots(1,1)
stackdat=[stackscale1,stackscale2,stackscale3]
#print stackdat
plt.hist(stackdat,bins=n_bins,normed=0,histtype='bar',stacked=False,color=["red","green","blue"])
#plt.hist(stackscale1,bins=n_bins,normed=0,histtype='bar')
#print stackscale1
#print stackscale2


#ax.hist(stackscale1,bins='fd',normed=0,histtype='bar')
#yikes=[]
#for element in stackscale1:
#	yikes.append(log(0.10125090724393926*expon.pdf(element,scale=0.00259501748635))-log(0.8987490927560609*pdf(element,0.9112770963751875,0.00259501748635)))

#ax.hist(yikes,bins='fd',normed=0,histtype='bar')

#print ((0.10125090724393926*expon.pdf(0.012164167,scale=0.00259501748635))/((0.10125090724393926*expon.pdf(0.012164167,scale=0.00259501748635))+0.8987490927560609*pdf(0.012164167,0.9112770963751875,0.00259501748635)))

x=linspace(0,0.03,1000)
y=0.397*expon.pdf(x,scale=0.002522)
#y=0.0680043311226034*expon.pdf(x,scale=0.00307142520134)
z=[]
w=[]
for index,element in enumerate(x):
	z.append(0.602*pdf(element,1.49,0.002522))
#	z.append(0.8986069666175454*pdf(element,1.0133732049359685,0.00307142520134))
##	w.append(0.0333887022598513*pdf(element,3.8947261339109507,0.00307142520134))
	w.append(y[index]+z[index])
#ax.plot(x,y,lw=3,alpha=0.8)
#ax.plot(x,z,lw=3,alpha=0.8)

#ax.plot(x,w,lw=3,alpha=0.8)
##w=[]
##for index,ele in enumerate(z):
##	w.append(ele/y[index])
##	if ele<=y[index]:
##		print index

#print w.index(max(w))
#print z.index(max(z))
#w=z/y
#print w
##ax.plot(x,w,lw=3,alpha=0.8)
#plt.axvline(0.0108675,color='red',lw=3,alpha=0.8)



#print 0.8987490927560609*pdf(0.0038675,0.9112770963751875,0.00259501748635)
#print 0.10125090724393926*expon.pdf(0.0038675,scale=0.00259501748635)
#print 2*iqr(stackscale1)/(pow(len(stackscale1),1/3))
#print ptp(stackscale1)
#f=open(sys.argv[3],'r+w')
#f.write(str(x))
#f.write('\n')
#f.write(str(y))
#f.write('\n')
#f.write(str(z))
#f.write('\n')
#f.write(str(stackscale1))
#f.close()
plt.show()
