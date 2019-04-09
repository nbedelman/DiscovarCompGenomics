#!/usr/bin/env python

import sys

DstatOut=sys.argv[1]
outFile=sys.argv[2]

eratoClade=['Htel','Hhsa','Hera','Hhim','HeraHhimHyb','Hsar','Hdem']
melpomeneClade=['Ldor','Hbes','Hpar','Hnum','Htim','Hcyd','Hmel','Hbur','Hhec','Hele']

d=open(DstatOut,"r")
o=open(outFile,"w")
for line in d:
    if "best" in line:
        atts=line.split()
        outGroup=atts[1]
        S3=atts[2]
        S2=atts[3]
        S1=atts[4]
        Dstat=atts[5]
        StdErr=atts[6]
        Zval=float(atts[7])
        if outGroup=='Etal':
            if (S1 in eratoClade) and (S2 in eratoClade) and (S3 in eratoClade):
                clade='erato'
            elif (S1 in melpomeneClade) and (S2 in melpomeneClade) and (S3 in melpomeneClade):
                clade='melpomene'
            else:
                clade='cross'
            if abs(Zval)>10:
                sig=3
            elif abs(Zval)>6:
                sig=2
            elif abs(Zval)>3:
                sig=1
            else:
                sig=0
            o.write('''%s,%s,%s,%s,%f,%s,%s\n''' % (outGroup,S3,S2,S1,Zval,clade,sig))
d.close()
o.close()
