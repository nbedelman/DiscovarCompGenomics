#!/bin/env python

import sys


pslFile=sys.argv[1]
# outBed=sys.argv[2]
minLength=int(sys.argv[2])

def pslToBed(pslFile, minLength):
    p=open(pslFile, "r")
    # o=open(outFile, "w")
    for line in p:
        pAtts=line.split("\t")
        alnLength=int(pAtts[0])
        if alnLength >= minLength:
            querScaf=pAtts[9]
            querLength=pAtts[10]
            querStart=int(pAtts[11])
            querEnd=int(pAtts[12])
            querStrand=pAtts[8][0]
            tarScaf=pAtts[13]
            tarLength=int(pAtts[14])
            tarStart=int(pAtts[15])
            tarEnd=int(pAtts[16])
            tarStrand=pAtts[8][1]
            if tarStrand == "+":
                relStart=querStart-tarStart
                relEnd=querEnd+(tarLength-tarEnd)
                col="0,0,255"
            elif tarStrand=="-":
                relStart=querStart-(tarLength-tarEnd)
                relEnd=querEnd+tarStart
                col="245,145,50"
            if relStart<0:
                relStart=0
            if relEnd>querLength:
                relEnd=querLength
            print('''%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s''' %
            (querScaf,relStart,relEnd,tarScaf,"0",tarStrand,querStart,querEnd,col))
    p.close()
    #o.close()

pslToBed(pslFile,minLength)
