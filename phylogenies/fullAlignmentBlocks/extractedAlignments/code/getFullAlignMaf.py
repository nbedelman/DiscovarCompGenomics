#!/bin/env python

import sys

mafFile=sys.argv[1]
allowedGaps=float(sys.argv[2])

maf = open(mafFile, "r")
fullAlign=open(mafFile[:-4]+"_fullAlign.maf", "w")
fullAlign.write("##maf version=1 scoring=N/A \n\n")
scafMap=open(mafFile[:-4]+"_fullAlign.scafMap.bed","w")

thisEntry='a\n'
csvEntry=''
goodSoFar=True
for line in maf:
    if line.startswith("a"):
        #print numSeqs
        if thisEntry != 'a\n':
            #print thisEntry
            fullAlign.write(thisEntry+"\n")
            scafMap.write(csvEntry)
        thisEntry='a\n'
        csvEntry=''
        goodSoFar=True
    elif line.startswith("s"):
        if goodSoFar==True:
            seq=line.split()[6]
            pctGaps=float(seq.count("-"))/len(seq)
            print(pctGaps)
            if pctGaps<=allowedGaps:
                thisEntry+=line
                if csvEntry=='':
                    csvEntry='''%s\t%s\t%s\t.\t%s\n''' % (line.split()[1],line.split()[2], int(line.split()[2])+int(line.split()[3]), line.split()[4])
            else:
                print("killed")
                goodSoFar=False
                thisEntry='a\n'
                csvEntry=''
fullAlign.close()
maf.close()
scafMap.close()
