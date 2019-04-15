#!/bin/env python

import sys
import os
import random
from Bio import SeqIO

faDir=sys.argv[1]
outBase=sys.argv[2]
numIntros=sys.argv[3]
numFiles=int(sys.argv[4])
numEntries=int(sys.argv[5])
minSize=int(sys.argv[6])
maxSize=int(sys.argv[7])
try:
    additionalOptions=sys.argv[8]
except IndexError:
    additionalOptions=''


def makeFastaDict(directory,minSize,maxSize):
    fastaDict={}
    record=0
    for subdir, dirs, files in os.walk(directory):
        for faFile in files:
            faLines=[]
            i=SeqIO.index(directory+"/"+faFile,"fasta")
            try:
                firstRecord=i[i.keys().next()]
            except StopIteration:
                pass
            length=len(firstRecord)
            if ((length>minSize) and (length<maxSize) and (len(i)==7)):
                for rec in i.keys():
                    faLine=i[rec].name+" "+str(i[rec].seq)
                    faLines.append(faLine)
                fastaDict[record]=(length,faLines)
                record+=1
    return fastaDict

def writeNexusFiles(fastaDict,numFiles,numEntries,outBase,numIntros,additionalOptions):
    for out in range(numFiles):
        fileName='''%s.s%i.nexus'''%(outBase,out)
        o=open(fileName,"w")
        phyloOutput='''%s.s%i_outDir'''%(outBase,out)
        o.write('#NEXUS\n\nBEGIN data;\n\tDimensions ntax=7 nchar=50000;\n\
        Format datatype=dna symbols="ACTG" missing=- gap=?;\n\tMatrix\n')
        region=0
        toUse=random.sample(xrange(len(fastaDict.keys())),numEntries)
        for entry in toUse:
            o.write('''[locus%i,%i]\n'''%(region,fastaDict[entry][0]))
            for line in fastaDict[entry][1]:
                o.write('''%s\n'''%(line))
            region+=1
        o.write('''\n\n;End;\n\nBEGIN PHYLONET;\n\nMCMC_SEQ (all) -mr %s -pl 32 %s -dir %s;\n\nEND;''' % (numIntros, additionalOptions,phyloOutput))
        o.close()

fastaDict=makeFastaDict(faDir,minSize,maxSize)
writeNexusFiles(fastaDict,numFiles,numEntries,outBase,numIntros,additionalOptions)
