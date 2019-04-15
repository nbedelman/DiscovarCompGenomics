#!/usr/bin/env python

import sys
from Bio import SeqIO
from operator import itemgetter

#This only works for my version of H. erato demaphoon. tsv format modelled after john davey's Hmel2.transitions.tsv

fastaFile=sys.argv[1]
outFile=sys.argv[2]

f=SeqIO.parse(open(fastaFile, "r"),"fasta")
out=open(outFile,"w")
seqStats=[]
for record in f:
    length=len(record)
    chrom=int(record.name.split("_")[1][3:])
    scaf=int(record.name.split("_")[2])
    name=record.name
    seqStats.append([chrom,scaf,length,name])

seqStats=sorted(seqStats,key=itemgetter(0,1))

chromPos=1
currentChrom=1
for seq in seqStats:
    chrom=seq[0]
    scaf=seq[1]
    length=seq[2]
    name=seq[3]
    if seq[0]==currentChrom:
        out.write('''%i\t%i\t%i\t%i\t%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%s\t%i\t%i\n''' %
        (chrom,chromPos,chromPos+length-1,0,"ok",0,name,1,length,"+",name,length,name,chromPos,chromPos+length-1))
        chromPos=chromPos+length
    else:
        chromPos=1
        currentChrom=chrom
        out.write('''%i\t%i\t%i\t%i\t%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%s\t%i\t%i\n''' %
        (chrom,chromPos,chromPos+length-1,0,"ok",0,name,1,length,"+",name,length,name,chromPos,chromPos+length-1))
        chromPos=chromPos+length

out.close()
