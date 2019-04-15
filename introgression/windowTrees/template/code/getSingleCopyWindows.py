#!/bin/env python

#In this version, also include a windows bed file and put the maf blocks into a file based on which window they fall into.
import sys
import os
from Bio import SeqIO

mafFile=sys.argv[1]
windowBed=sys.argv[2]
refGenome=sys.argv[3]

def ucscToBed(scaf,start,end,lenDict):
    scafLength=lenDict[scaf]
    newStart=scafLength-int(end)+1
    newEnd=scafLength-int(start)+1
    return (newStart,newEnd)

def makeQueryDict(queryGenome):
    g=SeqIO.parse(queryGenome,"fasta")
    gDict={}
    for record in g:
        gDict[record.id]=len(record)
    return gDict

queryDict=makeQueryDict(refGenome)

maf = open(mafFile, "r")
thisEntry='a\n'
numSeqs=0
for line in maf:
    if line.startswith("a"):
        if numSeqs == 4:
            #if it is a single copy, figure out which window it belongs to
            if strand=='-':
                newStart=ucscToBed(scaf,start,end,queryDict)[0]
                newEnd=ucscToBed(scaf,start,end,queryDict)[1]
                start=newStart
                end=newEnd
            windowCMD='''echo %s %s %s |sed -e 's/ [ ]*/\t/g'|bedtools intersect -a - -b %s -wb|awk '{print $7}' ''' % (scaf, start, end, windowBed)
            findWindow=os.popen(windowCMD).read().split()
            for window in set(findWindow):
                windowFile=open("mafs/"+window+".maf", "a+")
                windowFile.write(thisEntry+"\n")
                windowFile.close()
        thisEntry='a\n'
        numSeqs=0
        start=''
        end=''
        scaf=''
    elif line.startswith("s"):
        #if this is the first entry, record its features
        if thisEntry == 'a\n':
            atts=line.split()
            scaf=atts[1].split(".")[1]
            start=int(atts[2])
            end=start+int(atts[3])
            strand=atts[4]
        thisEntry+=line
        numSeqs+=1
maf.close()
