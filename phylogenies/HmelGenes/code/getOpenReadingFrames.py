#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio import Seq

aln=sys.argv[1]
out=sys.argv[2]
refGenome=sys.argv[3]


alignment=SeqIO.parse(aln,"fasta")

goodSeqs=[]
refSeqLength=0
for record in alignment:
    seqLen=len(record.seq.ungap("-"))
    remainder=seqLen%3
    if remainder==0:
        remainder=3
    record.seq=record.seq.ungap("-")+"N"*(3-remainder)
    readingFrame=Seq.translate(record.seq[:-3],table='Standard', stop_symbol='*', to_stop=False)
    if len(readingFrame.split("*"))==1:
        goodSeqs.append(record)
        if record.name==refGenome:
            refSeqLength=len(record.seq)
largeSeqs=[]
for record in goodSeqs:
    if float(len(record.seq))/refSeqLength >= .7:
        largeSeqs.append(record)
if len(largeSeqs)>=3:
    SeqIO.write(largeSeqs,out,"fasta")
