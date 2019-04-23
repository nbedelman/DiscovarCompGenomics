#!/usr/bin/env python

import sys
import os
sys.path.insert(0, '/Users/nbedelman/Documents/Mallet_lab/referenceScaffolding/ScaffoldingWithDiscovar/SWD-cli/SWD/comms')
scriptDirectory='/Users/nbedelman/Documents/Mallet_lab/DiscovarCompGenomics/scripts'
sys.path.insert(0,scriptDirectory)

from ContigClass import *
from getThinSegment import getThinSegment

bedFile=sys.argv[1]
mapFile=sys.argv[2]

def condenseBed(firstBedFile, outFile):
    con=Contig(firstBedFile)
    con.cullSegments()
    orderedSegs=con.orderSegs(con.getGoodSegments())
    combined=con.combineLists(orderedSegs,False)
    con.outputBed(combined,outFile)

def splitStrands(bedFile):
    b=open(bedFile,"r")
    p=open(bedFile+".plus.bed","w")
    n=open(bedFile+".minus.bed","w")
    for line in b:
        atts=line.split()
        strand=atts[5]
        if strand =='-':
            n.write(line)
        elif strand=='+':
            p.write(line)
    

def makeFirstBedFile(bedFile, outFile):
    b=open(bedFile,"r")
    bestLength=0
    bestAlign=''
    for line in b:
        atts=line.split()
        length=int(atts[7])-int(atts[6])
        if length > bestLength:
            bestAlign=line
            bestLength=length
    o=open(outFile,"w")
    b=open(bedFile,"r")
    o.write(bestAlign)
    for line in b:
        o.write(line)

# def getThinSegment(bedFile, outFile):
#     b=open(bedFile, "r")
#     o=open(outFile, "w")
#     downstreamPlus=[1000000000]
#     upstreamPlus=[0]
#     downstreamMinus=[1000000000]
#     upstreamMinus=[0]
#     for line in b:
#         atts=line.split()
#         strand=atts[5]
#         thickStart=int(atts[6])
#         thickEnd=int(atts[7])
#         thinStart=int(atts[1])
#         thinEnd=int(atts[2])
#         scaf=atts[0]
#         name=atts[3]
#         if strand=="+":
#             if thickStart<downstreamPlus[0]:
#                 downstreamPlus=[thickStart,thinStart,scaf,name]
#             if thickEnd>upstreamPlus[0]:
#                 upstreamPlus=[thickEnd,thinEnd]
#         elif strand=="-":
#             if thickStart<downstreamMinus[0]:
#                 downstreamMinus=[thickStart,thinStart,scaf,name]
#             if thickEnd>upstreamMinus[0]:
#                 upstreamMinus=[thickEnd,thinEnd]
#     try:
#         scaf=downstreamPlus[2]
#         name=downstreamPlus[3]
#         o.write('''%s\t%i\t%i\t%s\t0\t%s\t%i\t%i\t%s\n''' % (scaf,downstreamPlus[1],upstreamPlus[1],name,"+",downstreamPlus[1],upstreamPlus[1],"0,0,255"))
#     except IndexError:
#         scaf=downstreamMinus[2]
#         name=downstreamMinus[3]
#     try:
#         o.write('''%s\t%i\t%i\t%s\t0\t%s\t%i\t%i\t%s\n''' % (scaf,downstreamMinus[1],upstreamMinus[1],name,"-",downstreamMinus[1],upstreamMinus[1],"245,145,50"))
#     except IndexError:
#         pass
#     o.close()

splitStrands(bedFile)
try:
    makeFirstBedFile(bedFile+".plus.bed",bedFile+".plus.first.bed")
    condenseBed(bedFile+".plus.first.bed",bedFile+".plus.first.condense.bed")
    getThinSegment(bedFile+".plus.first.condense.bed",bedFile+".plus.first.condense.thin.bed")
except IndexError:
    pass

try:
    makeFirstBedFile(bedFile+".minus.bed",bedFile+".minus.first.bed")
    condenseBed(bedFile+".minus.first.bed",bedFile+".minus.first.condense.bed")
    getThinSegment(bedFile+".minus.first.condense.bed",bedFile+".minus.first.condense.thin.bed")
except IndexError:
    pass

os.system('''cat %s %s > %s ''' % (bedFile+".plus.first.condense.bed",bedFile+".minus.first.condense.bed",bedFile+".condensed.bed"))
os.system('''cat %s %s > %s ''' % (bedFile+".plus.first.condense.thin.bed",bedFile+".minus.first.condense.thin.bed",bedFile+".condensed.thin.bed"))

os.system('''%s/scafBedToChromBed.py %s %s %s''' % (scriptDirectory,bedFile+".condensed.bed",bedFile+".condensed.chroms.bed", mapFile))
os.system('''%s/scafBedToChromBed.py %s %s %s''' % (scriptDirectory,bedFile+".condensed.thin.bed",bedFile+".condensed.chroms.thin.bed", mapFile))
