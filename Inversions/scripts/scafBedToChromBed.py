#!/usr/bin/env python

import sys

bedFile=sys.argv[1]
outFile=sys.argv[2]
mapFile=sys.argv[3]


def makeConversionKey(mapFile):
    '''use a map file to create a key that converts scaffold coordinates to chromosomes'''
    convertDict={}
    m=open(mapFile, "r")
    #ignore header line
    m.readline()
    for line in m:
        atts=line.split()
        chrom=atts[0]
        chromStart=int(atts[1])
        chromEnd=int(atts[2])
        scafStart=int(atts[7])
        scafEnd=int(atts[8])
        scafName=atts[6]
        strand=atts[9]
        try:
            convertDict[scafName].append([chrom,chromStart,chromEnd,strand, scafStart,scafEnd])
        except KeyError:
            convertDict[scafName]=[[chrom,chromStart,chromEnd,strand, scafStart,scafEnd],]
    m.close()
    return convertDict


def scafBedToChromBed(bedFile,outFile,conversionKey):
    s=open(bedFile, "r")
    o=open(outFile, "w")
    scafDict={}
    startPiece=None
    endPiece=None
    for line in s:
        atts=line.split()
        scaf=atts[0]
        start=int(atts[1])
        end=int(atts[2])
        try:
            name=atts[3]
        except IndexError:
            name="."
        try:
            score=atts[4]
        except IndexError:
            score="."
        try:
            strand=atts[5]
        except IndexError:
            strand="+"
        try:
            thickStart=int(atts[6])
        except IndexError:
            thickStart=start
        except ValueError:
            thickStart=atts[6]
        try:
            thickEnd=int(atts[7])
        except IndexError:
            thickEnd=end
        except ValueError:
            thickEnd=atts[7]
        try:
            col=atts[8]
        except IndexError:
            col="."
        for piece in conversionKey[scaf]:
            if (start >= piece[-2]) and (start <= piece[-1]):
                startPiece=piece
                start=start-piece[-2]+1
            if (end >= piece[-2]) and (end <= piece[-1]):
                endPiece=piece
                end=end-endPiece[-2]+1
        if (startPiece==endPiece) and startPiece:
            chrom=startPiece[0]
            if startPiece[3]=="+":
                newStart=startPiece[1]+start-1
                newEnd=startPiece[1]+end-1
                if type(thickStart) is int:
                    newThickStart=startPiece[1]+thickStart-1
                    newThickEnd=startPiece[1]+thickEnd-1
                else:
                    newThickEnd=thickEnd
                    newThickStart=thickStart
            else:
                newEnd=startPiece[2]-start+1
                newStart=startPiece[2]-end+1
                if type(thickStart) is int:
                    newThickEnd=startPiece[2]-thickStart+1
                    newThickStart=startPiece[2]-thickEnd+1
                else:
                    newThickEnd=thickEnd
                    newThickStart=thickStart
                if strand=="+":
                    strand="-"
                else:
                    strand="+"
            o.write('''%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n''' % ("chr"+chrom,newStart,newEnd,name,score,strand,str(newThickStart),str(newThickEnd), col))
    o.close()
    s.close()

conversionKey=makeConversionKey(mapFile)
scafBedToChromBed(bedFile,outFile,conversionKey)
