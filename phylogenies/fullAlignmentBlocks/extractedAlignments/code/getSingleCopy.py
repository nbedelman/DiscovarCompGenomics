#!/bin/env python

import sys

mafFile=sys.argv[1]

maf = open(mafFile, "r")
singleCopy=open(mafFile[:-4]+"_singleCopy.maf", "w")
singleCopy.write("##maf version=1 scoring=N/A \n\n")
scafMap=open(mafFile[:-4]+"_scafMap.bed","w")

thisEntry='a\n'
csvEntry=''
numSeqs=0
specs=[]
for line in maf:
    if line.startswith("a"):
        #print numSeqs
        if ((numSeqs == len(set(specs))) and (numSeqs > 2)):
            #print thisEntry
            singleCopy.write(thisEntry+"\n")
            scafMap.write(csvEntry)
        thisEntry='a\n'
        csvEntry=''
        numSeqs=0
        specs=[]
    elif line.startswith("s"):
        thisEntry+=line
        specs.append(line.split()[1].split(".")[0])
        if csvEntry=='':
            csvEntry='''%s\t%s\t%s\t.\t%s\n''' % (line.split()[1],line.split()[2], int(line.split()[2])+int(line.split()[3]), line.split()[4])
        numSeqs+=1
singleCopy.close()
maf.close()
scafMap.close()
