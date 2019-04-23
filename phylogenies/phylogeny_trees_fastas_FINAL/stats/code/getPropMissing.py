#!/usr/bin/env python

import sys
from Bio import SeqIO

inFasta=sys.argv[1]

def getPropMissing(faFile):
    f=SeqIO.parse(faFile,"fasta")
    tot=0
    missing=0
    for record in f:
        tot+=len(record)
        thisMiss=record.seq.count("*")+record.seq.count("-")+record.seq.count("N")+record.seq.count("n")
        missing+=thisMiss
    print(float(missing)/tot)

getPropMissing(inFasta)
