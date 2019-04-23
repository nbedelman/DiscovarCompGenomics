#!/bin/env python

import sys

mafFile=sys.argv[1]
fullSeqNum=int(sys.argv[2])
outBase=sys.argv[3]

maf = open(mafFile, "r")
fullAlign=open(outBase+".fullAlign.maf", "w")
scafMap=open(outBase+".fullAlign_scafMap.bed","w")

thisEntry='a\n'
csvEntry=''
numSeqs=0
specs=[]
for line in maf:
    if line.startswith("a"):
        #print numSeqs
        if (numSeqs >= fullSeqNum):
            #print thisEntry
            fullAlign.write(thisEntry+"\n")
            scafMap.write(csvEntry)
        thisEntry='a\n'
        csvEntry=''
        numSeqs=0
        specs=[]
    elif line.startswith("s"):
        thisEntry+=line
        specs.append(line.split()[1].split(".")[0])
        if csvEntry=='':
            csvEntry='''%s\t%s\t%s\t.\t%s\n''' % (line.split()[1].split(".")[1],line.split()[2], int(line.split()[2])+int(line.split()[3]), line.split()[4])
        numSeqs+=1
fullAlign.close()
maf.close()
scafMap.close()
