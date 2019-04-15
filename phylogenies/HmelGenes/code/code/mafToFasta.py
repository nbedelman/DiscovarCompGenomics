#!/usr/bin/env python

import sys

mafFile=sys.argv[1]
outFile=sys.argv[2]

def mafToPhylip(mafFile, outFile):
    m=open(mafFile, "r")
    o=open(outFile, "w")
    seqDict={}
    seqsThisBlock=[]
    basesThisBlock=0
    overallSeqs=[]
    totalBases=0
    previousBases=""
    for line in m:
        if line.startswith("a"):
            for spec in overallSeqs:
                if spec not in seqsThisBlock:
                    seqDict[spec]+="?"*basesThisBlock
            totalBases+=basesThisBlock
            seqsThisBlock=[]
            basesThisBlock=0
        elif line.startswith("s"):
            atts=line.strip().split()
            spec=atts[1].split(".")[0]
            seq=atts[6]
            seqsThisBlock.append(spec)
            if basesThisBlock==0:
                basesThisBlock=len(seq)
            try:
                seqDict[spec]+=seq
            except KeyError:
                previousBases+="-"*totalBases
                seqDict[spec]=previousBases+seq
    for s in seqDict.keys():
        o.write('''>%s\n%s\n''' % (s,seqDict[s]))
    m.close()
    o.close()

mafToPhylip(mafFile,outFile)
