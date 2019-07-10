#!/bin/env python

import sys
from Bio import SeqIO

def filterGenome(fastaFile,cutoff, output):
    f=SeqIO.parse(open(fastaFile, "r"),"fasta")
    cutoff=int(cutoff)
    out=open(output,"a+")
    for record in f:
        if len(record)>=cutoff:
            SeqIO.write([record,],out,"fasta")
    f.close()
    out.close()

filterGenome(sys.argv[1],sys.argv[2],sys.argv[3])
