#!/usr/bin/env python

import sys

mafFile=sys.argv[1]
outFile=sys.argv[2]

def mafToPhylip(mafFile, outFile):
    m=open(mafFile, "r")
    o=open(outFile,"w")
    seqDict={}
    basesThisBlock=0
    for line in m:
        if line.startswith("a"):
            seqsThisBlock=[]
            basesThisBlock=0
        elif line.startswith("s"):
            atts=line.strip().split()
            spec=atts[1].split(".")[0]
            seq=atts[6].replace("-","?")
            seqsThisBlock.append(spec)
            if basesThisBlock==0:
                basesThisBlock=int(atts[3])
            seqDict[spec]=seq
    numSpecs=len(seqDict.keys())
    o.write('''\t%i\t%i\n\n''' % (numSpecs,basesThisBlock))
    for s in seqDict.keys():
        o.write('''^%s   %s\n''' % (s,seqDict[s]))
    m.close()
    o.close()

mafToPhylip(mafFile,outFile)
