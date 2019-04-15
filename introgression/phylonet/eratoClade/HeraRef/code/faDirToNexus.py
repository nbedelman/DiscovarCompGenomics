#!/bin/env python
#for window size 10KB; taking from analysis with sliding windows every 5KB


import sys
import os
from Bio import SeqIO

faDir=sys.argv[1]
output=sys.argv[2]
numIntros=sys.argv[3]

def faFileToNexEntry(faFile):
    try:
        num=int(faFile.split("_")[1].split(".")[0])
        if num%10 != 0 :
            return ""
    except IndexError:
        return ""
    faLines=[]
    i=SeqIO.parse(faFile,"fasta")
    try:
        firstRecord=i.next()
    except StopIteration:
        return ""
    length=len(firstRecord)
    if length>5000 and length < 15000:
        faLine=firstRecord.name+" "+str(firstRecord.seq)
        faLines.append(faLine)
        for record in i:
            faLine=record.name+" "+str(record.seq)
            faLines.append(faLine)
        return(length,faLines)
    else:
        return ""



def writeNexus(directory,output,numIntros):
    o=open(output,"w")
    phyloOutput=output[:-6]+"_outDir"
    o.write('#NEXUS\n\nBEGIN data;\n\tDimensions ntax=7 nchar=15000;\n\t\
    Format datatype=dna symbols="ACTG" missing=- gap=?;\n\tMatrix\n')
    region=0
    for subdir, dirs, files in os.walk(directory):
            for faFile in files:
                nexEntry=faFileToNexEntry(os.path.join(directory,faFile))
                if nexEntry != "":
                    o.write('''[locus%i,%i]\n'''%(region,nexEntry[0]))
                    for entry in nexEntry[1]:
                        o.write('''%s\n'''%(entry))
                region+=1
    o.write('''\n\n;End;\n\nBEGIN PHYLONET;\n\nMCMC_SEQ (all) -mr %s -pl 32 -dir %s;\n\nEND;''' % (numIntros, phyloOutput))
    o.close()



writeNexus(faDir,output,numIntros)
