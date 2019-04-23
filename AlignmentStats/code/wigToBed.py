#!/usr/bin/env python

#wigToBed.py

#Takes a wig file, and extracts regions of interest based on depth,length, and gap length
#Outputs a bed file with each high quality region

#Usage: wigToBed

import sys

def extractHighDepthBlocks(inputWig, outputBed, minCov, maxGap, minLength):
    f=open(inputWig, "r")
    beds=[]
    inBlock=False
    for line in f:
        if "fixed" in line:
            if inBlock:
                blockEnd=count-1
                beds.append([chrom, blockStart, blockEnd])
            inBlock=False
            inGap=False
            inTrial=False
            count=1
            chromPart=line.split()[1]
            chrom=chromPart.split("=")[1]
            trialStart=0
            blockStart=0
            gapStart=0
            blockEnd=0
        else:
            if not (inBlock or inGap or inTrial):
                #print '''in nothing, count = %s''' % (count)
                if int(line) >= minCov:
                    if minLength==1:
                        inBlock=True
                        blockStart=count-int(maxGap/2)
                        count+=1
                    else:
                        inTrial=True
                        trialStart=count
                        count+=1
                else:
                    count+=1
            elif inTrial:
                #print '''in trial, count = %s''' % (count)
                if int(line) >= minCov and ((count-trialStart)>=minLength):
                    inBlock=True
                    inTrial=False
                    blockStart=trialStart-int(maxGap/2)
                    count+=1
                elif int(line) >= minCov and ((count-trialStart)<=minLength):
                    count+=1
                else:
                    inTrial=False
                    count+=1
            elif inBlock and not inGap:
                #print '''in block, count = %s''' % (count)
                if int(line) >= minCov:
                    count+=1
                elif maxGap==0:
                    blockEnd=count-1
                    beds.append([chrom, blockStart, blockEnd])
                    inGap=False
                    inBlock=False
                    inTrial=False
                    count+=1
                else:
                    inGap=True
                    gapStart=count
                    count+=1
            elif inGap:
                if int(line) >= minCov:
                    inGap=False
                    inBlock=True
                    count+=1
                else:
                    if count-gapStart>=maxGap:
                        blockEnd=count-(maxGap/2)
                        beds.append([chrom, blockStart, blockEnd])
                        inGap=False
                        inBlock=False
                        inTrial=False
                        count+=1
                    else:
                        count+=1
    #if the block is still going at the end of the file, add it
    if inBlock:
        blockEnd=count
        beds.append([chrom, blockStart, blockEnd])
    #write the bed file
    out=open(outputBed, "w")
    for i in beds:
        out.write(str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\n")

inWig=sys.argv[1]
outBed=sys.argv[2]
minCov=int(sys.argv[3])
maxGap=int(sys.argv[4])
minLength=int(sys.argv[5])

extractHighDepthBlocks(inWig, outBed, minCov, maxGap, minLength)
