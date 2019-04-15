#!/usr/bin/env python

import sys

cnOut=sys.argv[1]
outFile=sys.argv[2]

def parseCmpNetOut(cnOut,outFile):
    c=open(cnOut,"r")
    o=open(outFile,"w")
    currentComp=""
    for line in c:
        if "Cmpnets" in line:
            atts=line.split()
            currentComp='''%s\t%s\t''' % (atts[1],atts[2])
        elif "distance" in line:
            score=line.split()[-1]
            o.write('''%s%s\n''' % (currentComp,score))
    c.close()
    o.close()

parseCmpNetOut(cnOut,outFile)
