#!/bin/env python

import sys
import os
import random
from Bio import SeqIO

faDir=sys.argv[1]
output=sys.argv[2]
numIntros=sys.argv[3]
totalAvailable=int(sys.argv[4])
totalDesired=int(sys.argv[5])
minSize=sys.argv[6]
try:
    additionalOptions=sys.argv[7]
except IndexError:
    additionalOptions=''

def faFileToNexEntry(faFile,seqNum,toUse, minSize):
    # num=faFile.split("_")[1].split(".")[0]
    faLines=[]
    i=SeqIO.parse(faFile,"fasta")
    try:
        firstRecord=i.next()
    except StopIteration:
        return ""
    length=len(firstRecord)
    if (length>int(minSize)):
        if seqNum in toUse:
            faLine=firstRecord.name+" "+str(firstRecord.seq)
            faLines.append(faLine)
            for record in i:
                faLine=record.name+" "+str(record.seq)
                faLines.append(faLine)
            return(length,faLines)
        else:
            return "skip"
    else:
        return ""



def writeNexus(directory,output,numIntros, totalAvailable,totalDesired, additionalOptions,minSize):
    o=open(output,"w")
    phyloOutput=output[:-6]+"_outDir"
    o.write('#NEXUS\n\nBEGIN data;\n\tDimensions ntax=7 nchar=50000;\n\t\
    Format datatype=dna symbols="ACTG" missing=- gap=?;\n\tMatrix\n')
    region=0
    seqNum=0
    toUse=random.sample(xrange(totalAvailable),totalDesired)
    for subdir, dirs, files in os.walk(directory):
            for faFile in files:
                nexEntry=faFileToNexEntry(os.path.join(directory,faFile), seqNum,toUse,minSize)
                if nexEntry=="skip":
                    seqNum+=1
                elif nexEntry != "":
                    seqNum+=1
                    o.write('''[locus%i,%i]\n'''%(region,nexEntry[0]))
                    for entry in nexEntry[1]:
                        o.write('''%s\n'''%(entry))
                    region+=1
    o.write('''\n\n;End;\n\nBEGIN PHYLONET;\n\nMCMC_SEQ (all) -mr %s -pl 32 %s -dir %s;\n\nEND;''' % (numIntros, additionalOptions,phyloOutput))
    o.close()


writeNexus(faDir,output,numIntros, totalAvailable,totalDesired, additionalOptions,minSize)
