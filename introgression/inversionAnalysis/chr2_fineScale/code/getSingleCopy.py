#!/bin/env python

import sys

mafFile=sys.argv[1]

maf = open(mafFile, "r")
singleCopy=open(mafFile[:-4]+"_singleCopy.maf", "w")
singleCopy.write("##maf version=1 scoring=N/A \n\n")
# scafMap=open(mafFile[:-4]+"scafMap.csv","w")

thisEntry='a\n'
# csvEntry=''
numSeqs=0
specs=[]
for line in maf:
    if line.startswith("a"):
        #print numSeqs
        if numSeqs == len(set(specs)):
            #print thisEntry
            singleCopy.write(thisEntry+"\n")
            # scafMap.write(csvEntry+"\n")
        thisEntry='a\n'
        # csvEntry=''
        numSeqs=0
        specs=[]
    elif line.startswith("s"):
        thisEntry+=line
        specs.append(line.split()[1].split(".")[0])
        # csvEntry+='''%s,%s,%s,%s,''' % (line.split()[1],line.split()[2], int(line.split()[2])+int(line.split()[3]), line.split()[4])
        numSeqs+=1
singleCopy.close()
maf.close()
# scafMap.close()
