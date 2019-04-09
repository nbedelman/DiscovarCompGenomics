#!/usr/bin/env python

import sys

snpFile=sys.argv[1]
genoFile=sys.argv[2]
outBase=sys.argv[3]
mapFile=sys.argv[4]


def makeConversionKey(mapFile):
    '''use a map file to create a key that converts scaffold coordinates to chromosomes'''
    convertDict={}
    m=open(mapFile, "r")
    #ignore header line
    m.readline()
    for line in m:
        atts=line.split()
        chrom=atts[0]
        scafStart=int(atts[1])
        scafEnd=int(atts[2])
        scafName=atts[6]
        strand=atts[9]
        convertDict[scafName]=[chrom,scafStart,scafEnd,strand]
    m.close()
    return convertDict


def scafSnpToChromSnp(snpFile,genoFile,outBase,conversionKey):
    s=open(snpFile, "r")
    g=open(genoFile, "r")
    snpOut=open(outBase+".snp", "w")
    genoOut=open(outBase+".geno", "w")
    scafDict={}
    for line in s:
        geno=g.readline()
        atts=line.split()
        scaf=atts[1]
        position=int(atts[3])
        chrom=conversionKey[scaf][0]
        if chrom != "0":
            if conversionKey[scaf][3]=="+":
                newPosition=conversionKey[scaf][1]+position-1
            else:
                newPosition=conversionKey[scaf][2]-position+1
            snpOut.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (atts[0],chrom,atts[2],newPosition,atts[4],atts[5]))
            genoOut.write(geno)
    g.close()
    s.close()
    snpOut.close()
    genoOut.close()

conversionKey=makeConversionKey(mapFile)
scafSnpToChromSnp(snpFile,genoFile,outBase,conversionKey)
