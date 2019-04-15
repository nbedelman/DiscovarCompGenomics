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

def pickRegion(windowList):
    return random.choice(windowList)

def pickRandomWindow(bedFile,outFile):
    out=open(outFile,"w")
    regionDict=makeRegionDict(bedFile)
    for window in regionDict.keys():
        out.write(pickRegion(regionDict[window]))
    out.close()

pickRandomWindow(inFile,outFile)
