#!/usr/bin/env python

#getThinSegment.py

import sys

fullBed=sys.argv[1]
outFile=sys.argv[2]

def getThinSegment(bedFile, outFile):
    b=open(bedFile, "r")
    o=open(outFile, "w")
    downstreamPlus=[1000000000]
    upstreamPlus=[0]
    downstreamMinus=[1000000000]
    upstreamMinus=[0]
    for line in b:
        atts=line.split()
        strand=atts[5]
        thickStart=int(atts[6])
        thickEnd=int(atts[7])
        thinStart=int(atts[1])
        thinEnd=int(atts[2])
        scaf=atts[0]
        name=atts[3]
        if strand=="+":
            if thickStart<downstreamPlus[0]:
                downstreamPlus=[thickStart,thinStart,scaf,name]
            if thickEnd>upstreamPlus[0]:
                upstreamPlus=[thickEnd,thinEnd]
        elif strand=="-":
            if thickStart<downstreamMinus[0]:
                downstreamMinus=[thickStart,thinStart,scaf,name]
            if thickEnd>upstreamMinus[0]:
                upstreamMinus=[thickEnd,thinEnd]
    try:
        scaf=downstreamPlus[2]
        name=downstreamPlus[3]
        o.write('''%s\t%i\t%i\t%s\t0\t%s\t%i\t%i\t%s\n''' % (scaf,downstreamPlus[1],upstreamPlus[1],name,"+",downstreamPlus[1],upstreamPlus[1],"0,0,255"))
    except IndexError:
        scaf=downstreamMinus[2]
        name=downstreamMinus[3]
    try:
        o.write('''%s\t%i\t%i\t%s\t0\t%s\t%i\t%i\t%s\n''' % (scaf,downstreamMinus[1],upstreamMinus[1],name,"-",downstreamMinus[1],upstreamMinus[1],"245,145,50"))
    except IndexError:
        pass
    o.close()

# getThinSegment(fullBed,outFile)
