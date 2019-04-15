#!/usr/bin/env python

#take a bed file with multiple entries per region, and pick one per region at random to use

import sys
import random

inFile=sys.argv[1]
outFile=sys.argv[2]



def makeRegionDict(bedFile):
    i=open(bedFile, "r")
    regionDict={}
    for line in i:
        atts=line.split()
        window=atts[3]
        try:
            regionDict[window].append(line)
        except KeyError:
            regionDict[window]=[line,]
    i.close()
    return regionDict

def pickRandomWindow(bedFile,outFile):
    out=open(outFile,"w")
    regionDict=makeRegionDict(bedFile)
    for window in regionDict.keys():
        num=0
        for entry in regionDict[window]:
            atts=entry.split()
            out.write('''%s\t%s\t%s\t%s_%i\t%s\n''' % (atts[0],atts[1],atts[2],atts[3],num,atts[4]))
            num+=1
    out.close()

pickRandomWindow(inFile,outFile)
